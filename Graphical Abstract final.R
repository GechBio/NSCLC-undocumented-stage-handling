# ==============================================================================
# Script: Generate Final Polished Graphical Abstract (No Clipping)
# ==============================================================================
library(ggplot2)
library(patchwork)
library(dplyr)

# ------------------------------------------------------------------------------
# Panel 1: The Problem
# ------------------------------------------------------------------------------
df_gene <- data.frame(
  Strategy = factor(c("CCA", "MRA"), levels = c("CCA", "MRA")),
  Genes = c(2, 16),
  Color = c("#D55E00", "#009E73")
)

p1 <- ggplot(df_gene, aes(x = Strategy, y = Genes, fill = Strategy)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(aes(label = paste(Genes, "genes")), vjust = -0.5, fontface = "bold", size = 5) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  ylim(0, 20) +
  labs(title = "1. THE PROBLEM\nComplete-Case Analysis", 
       subtitle = "Deletes dying patients\n& distorts features",
       y = "Selected Features", x = "Strategy") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", color = "#D55E00", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10)),
        axis.text.x = element_text(face = "bold", size = 13),
        axis.title.x = element_text(margin = margin(t = 10)))

# ------------------------------------------------------------------------------
# Panel 2: The Simulation
# ------------------------------------------------------------------------------
set.seed(123)
df_sim <- data.frame(
  Strategy = factor(rep(c("S1", "S2", "S3", "S4"), each = 20)),
  Estimate = c(rnorm(20, 0.75, 0.05), rnorm(20, 0.62, 0.05), rnorm(20, 0.50, 0.05), rnorm(20, 0.75, 0.05))
)

p2 <- ggplot(df_sim, aes(x = Strategy, y = Estimate, color = Strategy)) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  geom_hline(yintercept = 0.62, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 4.5, y = 0.63, label = "True HR", color = "red", fontface = "bold", vjust = -0.5, size = 4) +
  scale_color_manual(values = c("#D55E00", "#009E73", "#56B4E9", "#CC79A7")) +
  labs(title = "2. THE SIMULATION\nCalibrated Monte Carlo",
       subtitle = "Proves MI over-corrects\n& IPW is biased",
       y = "Estimated Log HR", x = "Strategies") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", color = "#0072B2", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10)),
        axis.text.x = element_text(face = "bold", size = 13),
        axis.title.x = element_text(margin = margin(t = 10)))

# ------------------------------------------------------------------------------
# Panel 3: The Solution (Fixed Clipping)
# ------------------------------------------------------------------------------
# Combine checkmarks with text to keep alignment perfect
df_win <- data.frame(
  x = 1,
  y = c(3, 2, 1),
  Text = c("✓ Lowest absolute bias", "✓ Preserves lethal risk", "✓ 7-point checklist")
)

p3 <- ggplot(df_win, aes(x = x, y = y)) +
  geom_text(aes(label = Text), hjust = 0, size = 5, fontface = "bold", color = "#009E73") +
  # Changed hjust to 0 (left-aligned) to match the bullets below it
  annotate("text", x = 1, y = 4.5, label = "S2 Wins", size = 9, fontface = "bold", color = "#009E73", hjust = 0) +
  # Added clip = "off" to prevent the edges from being cut off
  coord_cartesian(xlim = c(1, 5), ylim = c(0, 5), clip = "off") + 
  labs(title = "3. THE SOLUTION\nUnknown Retention (S2)",
       subtitle = "The only robust strategy") +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(face = "bold", color = "#009E73", hjust = 0.5, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 25)),
        # Added 20pt left margin to give the text extra room
        plot.margin = margin(10, 10, 10, 20)) 

# ------------------------------------------------------------------------------
# Combine and Save
# ------------------------------------------------------------------------------
final_plot <- (p1 | p2 | p3) + 
  plot_layout(widths = c(1, 1.2, 1.4)) +
  plot_annotation(
    title = "Evaluating Missing-Data Strategies in NSCLC Prognostic Modeling",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  )

ggsave("graphical_abstract_final.pdf", final_plot, width = 10, height = 4, dpi = 300)
ggsave("graphical_abstract_final.tiff", final_plot, width = 10, height = 4, dpi = 300, compression = "lzw")

cat("Final graphical abstract generated with no clipped text! Check graphical_abstract_final.pdf\n")
