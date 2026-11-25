# ==============================================================================
# Script 06b: NSCLC Missing Data Visualization - Independent-Replicate Monte Carlo Figure
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Generate Figure 6 (two-panel Bias and 95% CI Coverage plot) from
#           the independent-replicate Monte Carlo simulation results (Script 05b).
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load & Format Simulation Summary Data
#   MODULE 2: Panel Construction
#       - Block 2.1: Define Shared Palette & Theme Settings
#       - Block 2.2: Panel A - Bias with 95% Monte Carlo Intervals
#       - Block 2.3: Panel B - 95% Confidence Interval Coverage
#   MODULE 3: Assembly & Export
#       - Block 3.1: Combine Panels via Patchwork & Save Publication Figures
#
# Save as:  06b_plot_simulation_independent.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_sim_mc_independent.rds
# Output:   Figure6_Simulation_Independent.pdf, Figure6_Simulation_Independent.png
#
# Dependencies: ggplot2, dplyr, patchwork
#
# REVISION NOTES (Nov 2026):
#   - Refactored into standardized Module -> Block architecture matching Script 05.
#   - Added input file validation and factor-level normalization so both short
#     codes ("S0") and full labels ("S0: Full data (reference)") map cleanly.
#   - Applied colorblind-safe Okabe-Ito palette and unified legend placement at
#     the bottom-right/side to match manuscript Figure 6 styling.
#   - Exported both vector PDF and 300-DPI PNG files for journal submission.
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

# --- Block 1.2: Load & Format Simulation Summary Data -------------------------
load_and_format_summary <- function(rds_path = "nsclc_sim_mc_independent.rds") {
  if (!file.exists(rds_path)) {
    stop("Required input file not found: ", rds_path)
  }
  
  sim_obj <- readRDS(rds_path)
  tab <- if (is.list(sim_obj) && "summary" %in% names(sim_obj)) {
    sim_obj$summary
  } else {
    as.data.frame(sim_obj)
  }
  
  # Normalize strategy codes (extract leading S0, S1, S2, S3, S3b, S4 if long labels exist)
  tab$Strategy_Code <- sub(":.*$", "", trimws(as.character(tab$Strategy)))
  
  lv <- c("S0", "S1", "S2", "S3", "S3b", "S4")
  lb <- c(
    "S0\nFull data",
    "S1\nComplete-\ncase",
    "S2\nUnknown\nretained",
    "S3\nImputation,\nno outcome",
    "S3b\nImputation,\n+ outcome",
    "S4\nIPW"
  )
  
  tab$Strategy <- factor(tab$Strategy_Code, levels = lv, labels = lb)
  
  # Preserve standard mechanism order if MCAR/MAR/MNAR are present
  if (all(c("MCAR", "MAR", "MNAR") %in% unique(tab$Mechanism))) {
    tab$Mechanism <- factor(tab$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
  }
  
  return(tab)
}


# ==============================================================================
# MODULE 2: Panel Construction
# ==============================================================================

# --- Block 2.1: Define Shared Palette & Theme Settings ------------------------
get_shared_theme <- function() {
  theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0),
      axis.title.y     = element_text(size = 10.5, lineheight = 1.15),
      axis.text.x      = element_text(size = 9.5, colour = "black"),
      axis.text.y      = element_text(size = 9.5, colour = "black"),
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9.5),
      legend.background = element_blank()
    )
}

okabe_ito_palette <- c(
  "MCAR" = "#0072B2",  # Blue
  "MAR"  = "#009E73",  # Bluish green
  "MNAR" = "#CC79A7"   # Reddish purple
)

# --- Block 2.2: Panel A - Bias with 95% Monte Carlo Intervals -----------------
build_panel_bias <- function(tab, dodge_width = 0.6) {
  pos <- position_dodge(width = dodge_width)
  
  ggplot(tab, aes(x = Strategy, y = Bias_Adv, colour = Mechanism)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#D55E00", linewidth = 0.6) +
    geom_errorbar(
      aes(
        ymin = Bias_Adv - 1.96 * MCSE_bias,
        ymax = Bias_Adv + 1.96 * MCSE_bias
      ),
      width     = 0.22,
      linewidth = 0.65,
      position  = pos
    ) +
    geom_point(size = 2.6, position = pos) +
    scale_colour_manual(values = okabe_ito_palette, name = "Mechanism") +
    labs(
      title = "A. Bias",
      x     = NULL,
      y     = "Bias in Advanced-stage log HR\n(true value 0.62; bars = 95% Monte Carlo interval)"
    ) +
    get_shared_theme()
}

# --- Block 2.3: Panel B - 95% Confidence Interval Coverage --------------------
build_panel_coverage <- function(tab, dodge_width = 0.6, y_min = 0.70) {
  pos <- position_dodge(width = dodge_width)
  
  # Automatically expand lower y-limit if coverage drops below y_min
  min_cov <- min(tab$Coverage, na.rm = TRUE)
  lower_y <- if (is.finite(min_cov) && min_cov < y_min) max(0, floor(min_cov * 10) / 10) else y_min
  
  ggplot(tab, aes(x = Strategy, y = Coverage, colour = Mechanism)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "#D55E00", linewidth = 0.6) +
    geom_point(size = 2.8, position = pos) +
    scale_y_continuous(limits = c(lower_y, 1.00), breaks = seq(lower_y, 1.00, by = 0.05)) +
    scale_colour_manual(values = okabe_ito_palette, name = "Mechanism") +
    labs(
      title = "B. Coverage",
      x     = NULL,
      y     = "95% CI coverage"
    ) +
    get_shared_theme()
}


# ==============================================================================
# MODULE 3: Assembly & Export
# ==============================================================================

# --- Block 3.1: Combine Panels via Patchwork & Save Publication Figures -------
if (isTRUE(getOption("nsclc.plot_execute", TRUE))) {
  sim_summary <- load_and_format_summary("nsclc_sim_mc_independent.rds")
  
  p1 <- build_panel_bias(sim_summary)
  p2 <- build_panel_coverage(sim_summary)
  
  combined_fig6 <- (p1 + p2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
  
  ggsave(
    filename = "Figure6_Simulation_Independent.pdf",
    plot     = combined_fig6,
    width    = 11,
    height   = 4.8,
    device   = cairo_pdf
  )
  
  ggsave(
    filename = "Figure6_Simulation_Independent.png",
    plot     = combined_fig6,
    width    = 11,
    height   = 4.8,
    dpi      = 300
  )
  
  message("[Script 06b] Completed. Saved 'Figure6_Simulation_Independent.pdf' and '.png'.")
}
