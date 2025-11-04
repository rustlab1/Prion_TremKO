# Load libraries
library(limma)
library(ggplot2)

# Load data
counts <- read.csv("input_data/GSE281022_counts_with_annotation.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("input_data/GSE281022_metadata.csv", row.names = 1)

# Extract expression matrix and gene symbols
ex <- as.matrix(counts[, 5:ncol(counts)])  # assuming first 4 columns are annotation
genes <- counts$Gene.Symbol
rownames(ex) <- make.unique(genes)

# Match metadata order
ex <- ex[, rownames(meta)]

# Log2 transform if needed
if (max(ex, na.rm = TRUE) > 100) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

# Design and fit model
group <- factor(meta$Group, levels = c("WT", "Trem_KO"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(ex, design)
contrast <- makeContrasts(Trem_KO - WT, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))

# Get DE results
res <- topTable(fit2, adjust = "fdr", number = Inf)
res$gene <- rownames(res)

# Find top genes among significant ones (FDR < 0.05)
sig <- res[res$adj.P.Val < 0.05, ]
top_up <- sig[which.max(sig$logFC), "gene"]
top_down <- sig[which.min(sig$logFC), "gene"]

# Mark top genes for labeling and highlighting
res$label <- ifelse(res$gene %in% c(top_up, top_down), res$gene, "")
res$highlight <- ifelse(res$gene %in% c(top_up, top_down), "top", "other")

# ---- SAVE RESULTS ----
# Save full DE table
write.csv(res, "output.csv", row.names = FALSE)

# Save volcano plot as PNG
png("output.png", width = 1200, height = 900, res = 150)
ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = highlight), alpha = 0.6, size = ifelse(res$highlight == "top", 3, 1.5)) +
  scale_color_manual(values = c("top" = "red", "other" = "black"), guide = "none") +
  geom_text(
    data = subset(res, highlight == "top"),
    aes(label = label),
    vjust = -1.2, hjust = 0.5,
    color = "red", size = 4
  ) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 FDR") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  ylim(0, 7)
dev.off()
# -----------------------
