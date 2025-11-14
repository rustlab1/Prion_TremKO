#!/usr/bin/env python3
import argparse, hashlib, json, os, re, sys, uuid, tempfile, shutil, subprocess, csv
from pathlib import Path

INPUT_DIR = "input_data"
PAPER_NAME = "paper.pdf"

ALLOWED_ROOT_FILES = {
    "problem.json",
    "metadata.jsonl",
    "download.sh",
    "extracted_solution.csv",
    "extracted_solution.jpg",
    # plus exactly one of: plot.py or plot.R (handled separately)
}
# Disallowed names/patterns anywhere in the root
DISALLOWED_ROOT_DIRS = {
    ".git", "__pycache__", "node_modules", "venv", ".venv", "env", ".mamba", ".conda",
    ".Rproj.user", ".ipynb_checkpoints", ".pytest_cache", ".ruff_cache"
}
DISALLOWED_ROOT_FILES = {
    "Pipfile.lock", "poetry.lock", "package-lock.json", "yarn.lock",
    ".Rhistory", ".RData", ".Renviron", ".python-version", ".env",
    ".DS_Store", "requirements.txt",  # (requirements.txt implies multi-script env; disallow at root)
}
# Allowed *only* under input_data
ALLOWED_INPUT_EXTS = None  # any file allowed, we’ll verify via metadata.jsonl
REQUIRED_ROOT_FILES = ALLOWED_ROOT_FILES.copy()  # convenience

def md5sum(p: Path) -> str:
    h = hashlib.md5()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

def is_uuid_v4(name: str) -> bool:
    try:
        u = uuid.UUID(name, version=4)
        return str(u) == name
    except Exception:
        return False

def find_plot_file(root: Path):
    choices = [p for p in (root / "plot.py", root / "plot.R") if p.exists()]
    if not choices:
        return None
    if len(choices) > 1:
        return "BOTH"
    return choices[0]

def _read_text(path: Path):
    try:
        return path.read_text(errors="ignore")
    except Exception:
        return ""

def static_check_plot_io(plot_path: Path, strict: bool):
    if plot_path is None:
        return False, ["Missing plot.py/plot.R"], []
    if plot_path == "BOTH":
        return False, ["Both plot.py and plot.R present (choose exactly one)"], []

    txt = _read_text(plot_path)
    errors, warns = [], []

    # Must read extracted_solution.csv and only that as data input; must write extracted_solution.jpg
    in_ok  = "extracted_solution.csv" in txt
    out_ok = "extracted_solution.jpg" in txt
    if not in_ok:
        errors.append(f"{plot_path.name} must read 'extracted_solution.csv'.")
    if not out_ok:
        errors.append(f"{plot_path.name} must write 'extracted_solution.jpg'.")

    # Disallow hard-coded 'input_data/' (must assume same dir)
    if re.search(r'input_data[\\/]', txt):
        errors.append(f"{plot_path.name} must not reference 'input_data/...' (read from same directory).")

    # Heuristic: find other file-like refs
    suspicious = re.findall(
        r'["\']([^"\']+\.(csv|tsv|txt|xlsx|xls|json|rds|rda|png|jpg|jpeg|pdf))["\']',
        txt, flags=re.I
    )
    extra_refs = sorted(set([s[0] for s in suspicious
                             if s[0].lower() not in {"extracted_solution.csv", "extracted_solution.jpg"}]))
    if extra_refs:
        msg = f"{plot_path.name} references other files: {extra_refs}. It should read only the CSV and write only the JPG."
        if strict:
            errors.append(msg)
        else:
            warns.append(msg)

    return (len(errors) == 0), errors, warns

def read_metadata_jsonl(path: Path):
    errs, recs = [], []
    try:
        with path.open("r", encoding="utf-8") as f:
            for i, line in enumerate(f, 1):
                s = line.strip()
                if not s:
                    errs.append(f"metadata.jsonl line {i}: blank line not allowed.")
                    continue
                try:
                    obj = json.loads(s)
                except json.JSONDecodeError as e:
                    errs.append(f"metadata.jsonl line {i}: invalid JSON ({e}).")
                    continue
                if not isinstance(obj, dict):
                    errs.append(f"metadata.jsonl line {i}: must be a JSON object.")
                    continue
                if "name" not in obj or "md5sum" not in obj:
                    errs.append(f"metadata.jsonl line {i}: missing 'name' or 'md5sum'.")
                    continue
                if not isinstance(obj["name"], str) or not isinstance(obj["md5sum"], str):
                    errs.append(f"metadata.jsonl line {i}: 'name' and 'md5sum' must be strings.")
                    continue
                if not re.fullmatch(r"[0-9a-fA-F]{32}", obj["md5sum"]):
                    errs.append(f"metadata.jsonl line {i}: 'md5sum' must be 32 hex characters.")
                    continue
                recs.append(obj)
    except FileNotFoundError:
        errs.append("metadata.jsonl not found.")
    return recs, errs

def check_metadata_against_input(records, input_dir: Path):
    errs, warns = [], []
    actual_files = {p.name for p in input_dir.iterdir() if p.is_file() and not p.name.startswith(".")}
    meta_names = [r["name"] for r in records]
    meta_set = set(meta_names)

    missing = sorted(actual_files - meta_set)
    if missing:
        errs.append(f"metadata.jsonl is missing entries for: {missing}")

    extra = sorted(meta_set - actual_files)
    if extra:
        errs.append(f"metadata.jsonl lists files not found in {INPUT_DIR}/: {extra}")

    dups = sorted({n for n in meta_names if meta_names.count(n) > 1})
    if dups:
        errs.append(f"metadata.jsonl has duplicate entries: {dups}")

    for rec in records:
        name = rec["name"]
        path = input_dir / name
        if path.exists() and path.is_file():
            calc = md5sum(path)
            if calc.lower() != rec["md5sum"].lower():
                errs.append(f"MD5 mismatch for {name}: metadata={rec['md5sum']} actual={calc}")
    return errs, warns

def check_download_sh(script: Path):
    errs, warns = [], []
    if not script.exists():
        errs.append("download.sh is missing.")
        return errs, warns

    if not os.access(script, os.X_OK):
        warns.append("download.sh is not executable (chmod +x).")

    return errs, warns

def csv_is_valid_csv(path: Path):
    # Quick tests: extension, not XLSX/XLS magic, has commas, not only tabs
    if path.suffix.lower() != ".csv":
        return False, ["extracted_solution.csv must have .csv extension."]
    sig = path.read_bytes()[:8]
    if sig[:2] == b"PK":  # zip → likely xlsx
        return False, ["extracted_solution.csv appears to be an XLSX (zip) file."]
    if sig.startswith(b"\xD0\xCF\x11\xE0"):  # old xls OLE
        return False, ["extracted_solution.csv appears to be an Excel .xls file."]
    # Parse a few rows with the csv module (comma delimiter)
    try:
        with path.open("r", encoding="utf-8", newline="") as f:
            sniffer = csv.Sniffer()
            sample = f.read(2048)
            f.seek(0)
            dialect = None
            try:
                dialect = sniffer.sniff(sample)
            except Exception:
                pass
            if dialect and dialect.delimiter == "\t":
                return False, ["extracted_solution.csv looks like a TSV (tab-delimited). Use commas."]
            f.seek(0)
            reader = csv.reader(f)  # default comma
            # Read a few rows to ensure it parses
            for _ in range(5):
                next(reader)
    except StopIteration:
        # Empty CSV is still syntactically valid; allow
        return True, []
    except UnicodeDecodeError:
        return False, ["extracted_solution.csv is not UTF-8 text."]
    except Exception as e:
        return False, [f"extracted_solution.csv could not be parsed as CSV: {e}"]
    # Also basic comma presence check
    txt = path.read_text(errors="ignore")
    if ("," not in txt) and ("\n" in txt.strip()):
        return False, ["extracted_solution.csv has no commas—did you save as TSV?"]
    return True, []

def enforce_no_extras(root: Path, plot_path: Path):
    """Ensure there are no other files/folders at root besides the allowlist + plot script + input_data."""
    errs = []

    allowed = set(ALLOWED_ROOT_FILES)
    if plot_path and plot_path != "BOTH":
        allowed.add(plot_path.name)
    allowed.add(INPUT_DIR)

    root_items = [p.name for p in root.iterdir()]

    # Disallowed known dirs/files
    for dis in DISALLOWED_ROOT_DIRS:
        if dis in root_items:
            errs.append(f"Disallowed directory at root: {dis}")
    for dis in DISALLOWED_ROOT_FILES:
        if dis in root_items:
            errs.append(f"Disallowed file at root: {dis}")

    # Any unexpected items?
    unexpected = sorted([n for n in root_items if n not in allowed and not n.startswith(".")])
    if unexpected:
        errs.append(f"Unexpected items in root (only one script allowed and no extra files): {unexpected}")

    # Enforce “expert solution entirely in 1 script”: only plot.{py,R} as a script at root
    other_scripts = [n for n in root_items
                     if re.search(r'\.(py|r|ipynb|sh)$', n, flags=re.I)
                     and n.lower() not in {"download.sh"}
                     and (plot_path is None or n != plot_path.name)]
    if other_scripts:
        errs.append(f"Multiple scripts found at root; keep your expert solution in a single plot script. Extra scripts: {other_scripts}")

    # No subfolders except input_data
    other_dirs = [n for n in root_items if (root / n).is_dir() and n != INPUT_DIR]
    if other_dirs:
        errs.append(f"Extra folders at root (only '{INPUT_DIR}' is allowed): {other_dirs}")

    return errs

def static_check_scripts_no_inputdir_reads(root: Path):
    """Check all scripts at root (plot + download.sh) do not read from 'input_data/...'
       except download.sh which *should* target input_data writes; we only block plot.* here."""
    errs = []
    plot_py = root / "plot.py"
    plot_r  = root / "plot.R"
    for p in [plot_py, plot_r]:
        if p.exists():
            txt = _read_text(p)
            if re.search(r'input_data[\\/]', txt):
                errs.append(f"{p.name} must not read from 'input_data/...'. Use files in the same directory as the script.")
            # Also discourage absolute paths
            if re.search(r'["\']/(?!dev/null)[^"\']+', txt):
                errs.append(f"{p.name} appears to use absolute file paths; use local filenames only.")
    return errs

def try_run_plot(root: Path, plot_path: Path):
    """Run plot.{py,R} in a temp dir with only the plot file and extracted_solution.csv present.
       Success = process exits 0 and produces extracted_solution.jpg (and nothing else created)."""
    if plot_path is None or plot_path == "BOTH":
        return False, ["Cannot run plot script (missing or both present)."], []

    csv_src = root / "extracted_solution.csv"
    if not csv_src.exists():
        return False, ["extracted_solution.csv missing; cannot run plot test."], []

    errors, warns = [], []
    with tempfile.TemporaryDirectory() as tmpd:
        tmp = Path(tmpd)
        # Copy only these two files
        dest_plot = tmp / plot_path.name
        dest_csv  = tmp / "extracted_solution.csv"
        shutil.copy2(plot_path, dest_plot)
        shutil.copy2(csv_src, dest_csv)

        before = set(os.listdir(tmp))
        try:
            if plot_path.suffix.lower() == ".py":
                cmd = ["python3", dest_plot.name]
            else:
                # Prefer Rscript; fall back to R -f
                if shutil.which("Rscript"):
                    cmd = ["Rscript", dest_plot.name]
                else:
                    cmd = ["R", "-f", dest_plot.name]
            proc = subprocess.run(cmd, cwd=tmp, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
            if proc.returncode != 0:
                errors.append(f"Running {plot_path.name} failed (exit {proc.returncode}). Stderr:\n{proc.stderr.decode(errors='ignore')[:500]}")
        except FileNotFoundError as e:
            warns.append(f"Runtime tool not available to run {plot_path.name}: {e}. (Static checks still applied.)")
            return True, [], warns  # don’t fail if the runtime (R/Python) isn’t present
        except subprocess.TimeoutExpired:
            errors.append(f"Running {plot_path.name} timed out after 120s.")
        after = set(os.listdir(tmp))

        # Check outputs
        if "extracted_solution.jpg" not in after:
            errors.append(f"{plot_path.name} did not create 'extracted_solution.jpg' in the same directory.")
        # Ensure no unexpected outputs (besides common harmless __pycache__)
        created = sorted(list(after - before))
        extras = [x for x in created if x not in {"extracted_solution.jpg", "__pycache__"}]
        if extras:
            errors.append(f"{plot_path.name} created unexpected outputs: {extras} (should only write extracted_solution.jpg).")

    return (len(errors) == 0), errors, warns

def main():
    ap = argparse.ArgumentParser(description="Validate submission packaging & rules.")
    ap.add_argument("--strict", action="store_true", help="Fail if plot script references any extra files.")
    ap.add_argument("--run-plot", action="store_true", help="Attempt to execute plot.{py,R} in a clean temp dir.")
    ap.add_argument("submission_path", help="Path to the <UUIDv4> folder")
    args = ap.parse_args()

    root = Path(args.submission_path).resolve()
    report, failures = [], 0

    # 1) Folder name UUIDv4 good
    name_ok = is_uuid_v4(root.name)
    report.append(("Folder name is UUIDv4", name_ok, [] if name_ok else [f"Folder name '{root.name}' is not a valid UUIDv4."]))
    failures += 0 if name_ok else 1

    # 2) Required structure good
    input_dir = root / INPUT_DIR
    input_ok = input_dir.exists() and input_dir.is_dir()
    report.append((f"'{INPUT_DIR}/' subfolder exists", input_ok, [] if input_ok else [f"Missing {INPUT_DIR}/ subfolder."]))
    failures += 0 if input_ok else 1

    # 3) Required root files present good
    missing = [f for f in REQUIRED_ROOT_FILES if not (root / f).exists()]
    req_ok = (len(missing) == 0)
    report.append(("Required root files present", req_ok, [] if req_ok else [f"Missing required root files: {missing}"]))
    failures += 0 if req_ok else 1

    # 4) Exactly one plot.{py,R} good + static IO rules (+ no input_data/ path usage) good
    plot = find_plot_file(root)
    ok_plot, plot_errs, plot_warns = static_check_plot_io(plot, args.strict)
    report.append(("Plot script validity (IO rules, local paths)", ok_plot, plot_errs + [f"Warning: {w}" for w in plot_warns]))
    failures += 0 if ok_plot else 1
    script_path_errs = static_check_scripts_no_inputdir_reads(root)
    if script_path_errs:
        report.append(("Scripts do not read via 'input_data/...'", False, script_path_errs))
        failures += 1
    else:
        report.append(("Scripts do not read via 'input_data/...'", True, []))

    # 5) input_data contains paper.pdf (supplementary PDFs optional) good
    if input_ok:
        paper_ok = (input_dir / PAPER_NAME).is_file()
        report.append((f"{INPUT_DIR}/{PAPER_NAME} present", paper_ok, [] if paper_ok else [f"{PAPER_NAME} not found in {INPUT_DIR}/."]))
        failures += 0 if paper_ok else 1

    # 6) metadata.jsonl vs input_data + MD5s good
    metadata = root / "metadata.jsonl"
    if metadata.exists() and input_ok:
        recs, rec_errs = read_metadata_jsonl(metadata)
        ok_meta = (len(rec_errs) == 0)
        report.append(("metadata.jsonl structure valid", ok_meta, rec_errs))
        failures += 0 if ok_meta else 1
        if ok_meta:
            m_errs, m_warns = check_metadata_against_input(recs, input_dir)
            ok2 = (len(m_errs) == 0)
            report.append(("metadata.jsonl matches files & MD5s", ok2, m_errs + [f"Warning: {w}" for w in m_warns]))
            failures += 0 if ok2 else 1

    # 7) download.sh basic checks (writes confined to input_data/)
    dl_errs, dl_warns = check_download_sh(root / "download.sh")
    dl_ok = (len(dl_errs) == 0)
    report.append(("download.sh basic checks", dl_ok, dl_errs + [f"Warning: {w}" for w in dl_warns]))
    failures += 0 if dl_ok else 1

    # 8) problem.json is valid JSON
    pj = root / "problem.json"
    if pj.exists():
        try:
            _ = json.loads(pj.read_text())
            report.append(("problem.json is valid JSON", True, []))
        except Exception as e:
            report.append(("problem.json is valid JSON", False, [f"Invalid JSON: {e}"]))
            failures += 1

    # 9) extracted_solution.csv is a real CSV (not TSV/XLSX/etc.) good
    csv_path = root / "extracted_solution.csv"
    if csv_path.exists():
        ok_csv, csv_errs = csv_is_valid_csv(csv_path)
        report.append(("extracted_solution.csv is valid CSV", ok_csv, csv_errs))
        failures += 0 if ok_csv else 1

    # 10) No extra files/folders in root; expert solution is in exactly one plot script good
    extra_errs = enforce_no_extras(root, plot)
    extras_ok = (len(extra_errs) == 0)
    report.append(("No extra files/folders in root; single-script solution", extras_ok, extra_errs))
    failures += 0 if extras_ok else 1

    # Print report
    print("\nSUBMISSION CHECK REPORT\n" + "="*25)
    for title, ok, messages in report:
        status = "PASS" if ok else "FAIL"
        print(f"[{status}] {title}")
        for msg in messages:
            print(f"  - {msg}")
    print("="*25)
    if failures:
        print(f"Result: FAIL ({failures} check{'s' if failures!=1 else ''} failed)")
        sys.exit(1)
    else:
        print("Result: PASS")
        sys.exit(0)

if __name__ == "__main__":
    main()
