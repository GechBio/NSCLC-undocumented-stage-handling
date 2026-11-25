# ==============================================================================
# Script 03b: NSCLC Missing Data Visualization - 100-Split LASSO Selection Frequency (Fig 3)
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Generate Figure 3 (two-panel horizontal bar chart of LASSO gene
#           selection frequencies across 100 repeated stratified 70/30 splits
#           for Strategy A vs. Strategy B in the expression cohort).
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Repeated-Split Results (from Script 12)
#   MODULE 2: Data Preparation
#       - Block 2.1: Format & Order Candidate Gene Selection Frequencies
#   MODULE 3: Panel Construction & Export
#       - Block 3.1: Build Strategy A & Strategy B Horizontal Bar Panels
#       - Block 3.2: Combine via Patchwork & Save Figure 3 (PDF & 300-DPI PNG)
#
# Save as:  03b_generate_selection_frequency_plot.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_repeated_splits.rds
# Output:   Figure3_Selection_Frequency_100Splits.pdf (.png)
#
# Dependencies: ggplot2, dplyr, patchwork
#
# REVISION NOTES (Nov 2026):
#   - Replaced single-split dot plot (moved to Table 2) with the 100-split
#     selection frequency bar chart matching Figure 3 of the revised manuscript.
#   - Ordered candidate genes by combined selection frequency across strategies
#     and added vertical dashed reference lines at frequency = 0.50.
# ==============================================================================


# ==============================================================================
# MODULE 1: Environment & Dependency Management
# ==============================================================================

# --- Block 1.1: Load Required Packages ----------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# --- Block 1.2: Load Repeated-Split Results (from Script 12) ------------------
load_repeated_splits <- function(rds_path = "nsclc_repeated_splits.rds") {
  if (!file.exists(rds_path)) {
    stop("Run Script 12 first to generate ", rds_path)
  }
  readRDS(rds_path)
}


# ==============================================================================
# MODULE 2: Data Preparation
# ==============================================================================

# --- Block 2.1: Format & Order Candidate Gene Selection Frequencies -----------
prepare_selection_freq_data <- function(rep_out) {
  df_A <- data.frame(
    Gene     = names(rep_out$A$freq),
    Freq_A   = as.numeric(rep_out$A$freq),
    stringsAsFactors = FALSE
  )
  df_B <- data.frame(
    Gene     = names(rep_out$B$freq),
    Freq_B   = as.numeric(rep_out$B$freq),
    stringsAsFactors = FALSE
  )
  
  merged <- full_join(df_A, df_B, by = "Gene") %>%
    mutate(
      Freq_A     = ifelse(is.na(Freq_A), 0, Freq_A),
      Freq_B     = ifelse(is.na(Freq_B), 0, Freq_B),
      Mean_Freq  = (Freq_A + Freq_B) / 2
    ) %>%
    arrange(Mean_Freq, Freq_B, Gene) %>%
    mutate(Gene = factor(Gene, levels = Gene))
  
  return(merged)
}


# ==============================================================================
# MODULE 3: Panel Construction & Export
# ==============================================================================

# --- Block 3.1: Build Strategy A & Strategy B Horizontal Bar Panels -----------
build_freq_panel <- function(df, freq_col, panel_title, bar_color, show_y_labels = TRUE) {
  ggplot(df, aes(x = .data[[freq_col]], y = Gene)) +
    geom_col(fill = bar_color, width = 0.72) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "gray45", linewidth = 0.7) +
    scale_x_continuous(
      limits = c(0, 1.0),
      breaks = seq(0, 1.0, by = 0.2),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = panel_title,
      x     = "Selection frequency (100 splits)",
      y     = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", size = 11.5, hjust = 0.5),
      axis.title.x = element_text(size = 10.5, margin = margin(t = 8)),
      axis.text.x  = element_text(size = 9.5, color = "black"),
      axis.text.y  = if (show_y_labels) {
        element_text(face = "italic", size = 9.5, color = "black")
      } else {
        element_blank()
      }
    )
}

# --- Block 3.2: Combine via Patchwork & Save Figure 3 (PDF & 300-DPI PNG) -----
if (isTRUE(getOption("nsclc.fig03b_execute", TRUE))) {
  rep_out <- load_repeated_splits("nsclc_repeated_splits.rds")
  freq_df <- prepare_selection_freq_data(rep_out)
  
  p_left <- build_freq_panel(
    df            = freq_df,
    freq_col      = "Freq_A",
    panel_title   = "Strategy A (complete-case)",
    bar_color     = "#D55E00",
    show_y_labels = TRUE
  )
  
  p_right <- build_freq_panel(
    df            = freq_df,
    freq_col      = "Freq_B",
    panel_title   = "Strategy B (Unknown retained)",
    bar_color     = "#0072B2",
    show_y_labels = FALSE
  )
  
  fig3 <- (p_left + p_right) +
    plot_annotation(
      title = "LASSO selection frequency across 100 repeated 70/30 splits (expression)",
      theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0))
    )
  
  ggsave("Figure3_Selection_Frequency_100Splits.pdf", fig3, width = 9.5, height = 7.2, device = cairo_pdf)
  ggsave("Figure3_Selection_Frequency_100Splits.png", fig3, width = 9.5, height = 7.2, dpi = 300)
  
  message("[Script 03b] Completed. Saved 'Figure3_Selection_Frequency_100Splits.pdf' and '.png'.")
}
