library(ggplot2)

# Load the summary table (replace with your file)
deg_summary <- read.csv("input_data/data.csv", check.names = FALSE)

# Ensure correct order
deg_summary$Comparison <- factor(deg_summary$Comparison, levels = deg_summary$Comparison)

# Plot (identical style)
p <- ggplot(deg_summary, aes(x = Comparison, y = Pval005, fill = Comparison)) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = Pval005), vjust = -0.3, size = 5, color = "black") +
  geom_text(aes(label = FDR005, y = 10), vjust = 1.5, size = 5, color = "red") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "Comparison",
    y = "Number of DEGs p value < 0.005"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1, size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(size = 13)
  ) +
  scale_fill_manual(values = c("#619CFF", "#F8766D", "#D39200", "#FFD700")) +
  theme_bw()

# Display
p

# Save
ggsave("output_plot.png", plot = p, width = 8, height = 8)
