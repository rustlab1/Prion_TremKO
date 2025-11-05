# Load libraries
library(limma)
library(ggplot2)
library(patchwork)

# Load data
counts <- read.csv("input_data/GSE281022_annotated_expression_matrix_bar1_4.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("input_data/GSE281022_TremKO_prion_and_TremKO_NBH_metadata_bar1-4.csv", row.names = 1)
meta$Group <- gsub("\\.", "_", meta$Group)  # Normalize group names

# Extract numeric expression matrix and symbols
ex_raw <- as.matrix(sapply(counts[, grep("^GSM", colnames(counts))], as.numeric))
gene_symbols <- as.character(counts$Gene.symbol)

# Filter low-expression probes: RMA > 4 in at least 2 samples
keep <- rowSums(ex_raw > 4) >= 2
ex_filtered <- ex_raw[keep, ]
gene_symbols <- gene_symbols[keep]

# Collapse to gene level using limma::avereps
ex_gene <- avereps(ex_filtered, ID = gene_symbols)

# Initialize plot list
plots <- list()

### ---- Analysis block (reusable structure) ----
run_contrast_plot <- function(ex, meta, group_levels, contrast_label, plot_title) {
  meta_sub <- meta[meta$Group %in% group_levels, , drop = FALSE]
  ex_sub <- ex[, rownames(meta_sub)]
  
  group <- factor(meta_sub$Group, levels = group_levels)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(ex_sub, design)
  contrast <- makeContrasts(contrasts = paste(group_levels[2], "-", group_levels[1]), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  res <- topTable(fit2, adjust = "fdr", number = Inf)
  
  df <- data.frame(
    Category = c("P.Value < 0.005", "adj.P.Val < 0.05"),
    Count = c(sum(res$P.Value < 0.005), sum(res$adj.P.Val < 0.05))
  )
  
  p <- ggplot(df, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5, size = 4) +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10, face = "bold")) +
    labs(title = plot_title, x = "", y = "Significant genes")
  
  return(p)
}

# ---- Run all four comparisons ----
plots[[1]] <- run_contrast_plot(ex_gene, meta, c("TREM_NBH", "TREM_KO_ME7"), "TREM_KO_ME7 - TREM_NBH", "TREM_KO_ME7 vs TREM_NBH")
plots[[2]] <- run_contrast_plot(ex_gene, meta, c("WT_ME7", "TREM_KO_ME7"), "TREM_KO_ME7 - WT_ME7", "TREM_KO_ME7 vs WT_ME7")
plots[[3]] <- run_contrast_plot(ex_gene, meta, c("WT_NBH", "TREM_NBH"), "TREM_NBH - WT_NBH", "TREM_NBH vs WT_NBH")
plots[[4]] <- run_contrast_plot(ex_gene, meta, c("WT_NBH", "WT_ME7"), "WT_ME7 - WT_NBH", "WT_ME7 vs WT_NBH")

# ---- Combine all plots into one display ----
(plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])





# Re-run the contrast part only to get the DEG result tables
get_degs <- function(ex, meta, group_levels) {
  meta_sub <- meta[meta$Group %in% group_levels, , drop=FALSE]
  ex_sub <- ex[, rownames(meta_sub)]
  group <- factor(meta_sub$Group, levels = group_levels)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(ex_sub, design)
  contrast <- makeContrasts(contrasts = paste(group_levels[2], "-", group_levels[1]), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  topTable(fit2, adjust = "fdr", number = Inf)
}

# Store DEG tables
res1 <- get_degs(ex_gene, meta, c("TREM_NBH", "TREM_KO_ME7"))
res2 <- get_degs(ex_gene, meta, c("WT_ME7", "TREM_KO_ME7"))
res3 <- get_degs(ex_gene, meta, c("WT_NBH", "TREM_NBH"))
res4 <- get_degs(ex_gene, meta, c("WT_NBH", "WT_ME7"))

# Build DEG summary table
deg_summary <- data.frame(
  Comparison = c(
    "Trem2-/- prion vs Trem2-/- NBH",
    "Trem2-/- prion vs WT prion",
    "Trem2-/- NBH vs WT NBH",
    "WT prion vs WT NBH"
  ),
  Pval005 = c(
    sum(res1$P.Value < 0.005),
    sum(res2$P.Value < 0.005),
    sum(res3$P.Value < 0.005),
    sum(res4$P.Value < 0.005)
  ),
  FDR005 = c(
    sum(res1$adj.P.Val < 0.05),
    sum(res2$adj.P.Val < 0.05),
    sum(res3$adj.P.Val < 0.05),
    sum(res4$adj.P.Val < 0.05)
  )
)

deg_summary$Comparison <- factor(deg_summary$Comparison, levels = deg_summary$Comparison)


# Plot paper-style DEG barplot
aa <- ggplot(deg_summary, aes(x = Comparison, y = Pval005, fill = Comparison)) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = Pval005), vjust = -0.3, size = 5, color = "black") +
  geom_text(aes(label = FDR005, y = 10), vjust = 1.5, size = 5, color = "red") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "Comparison",
    y = "Number of DEGs (p value < 0.005)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1, size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(size = 13)
  ) +
  scale_fill_manual(values = c("#619CFF", "#F8766D", "#D39200", "#FFD700"))+
  theme_bw()
aa

# ---- SAVE DEG TABLES ----
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add each DEG result as a separate sheet
addWorksheet(wb, "TREM_KO_ME7 vs TREM_NBH")
writeData(wb, sheet = 1, x = res1, rowNames = TRUE)

addWorksheet(wb, "TREM_KO_ME7 vs WT_ME7")
writeData(wb, sheet = 2, x = res2, rowNames = TRUE)

addWorksheet(wb, "TREM_NBH vs WT_NBH")
writeData(wb, sheet = 3, x = res3, rowNames = TRUE)

addWorksheet(wb, "WT_ME7 vs WT_NBH")
writeData(wb, sheet = 4, x = res4, rowNames = TRUE)

# Save to file
saveWorkbook(wb, file = "output.xlsx", overwrite = TRUE)

# ---- SAVE Plot TABLES ----
ggsave("output.png", plot = aa, width = 8, height = 8)

