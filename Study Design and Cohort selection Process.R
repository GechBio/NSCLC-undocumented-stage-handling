# ==============================================================================
# Script 01b: NSCLC Missing Data Visualization - Study Design Flowchart (Figure 1)
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Generate Figure 1 (dual-arm study design flowchart) in R (ggplot2),
#           matching the revised independent-replicate simulation design
#           (S1, S2, S3, S3b, S4) and updated empirical/simulation metrics.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#   MODULE 2: Flowchart Geometry Definition
#       - Block 2.1: Define Nodes (Boxes, Text Labels, Fills, and Borders)
#       - Block 2.2: Define Directed Connectors (Arrows)
#   MODULE 3: Rendering & Export
#       - Block 3.1: Build ggplot2 Flowchart & Save Publication Files (PDF/PNG)
#
# Save as:  01b_generate_study_flow_diagram.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    None (self-contained schematic specification)
# Output:   Figure1_Final.pdf, Figure1_Final.png
#
# Dependencies: ggplot2
#
# REVISION NOTES (Nov 2026):
#   - Updated boxes and labels to match final manuscript Figure 1:
#       * Replaced "Pure-MAR Control" with "S3b: Imputation (+ outcome)"
#       * Added "(survival times regenerated in each replicate)" to missingness box
#       * Updated empirical & simulation output labels to match Table 1 & Table 4
#   - Refactored into standardized Module -> Block architecture.
# ==============================================================================


# ==============================================================================
# MODULE 1: Environment & Dependency Management
# ==============================================================================

# --- Block 1.1: Load Required Packages ----------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
})


# ==============================================================================
# MODULE 2: Flowchart Geometry Definition
# ==============================================================================

# --- Block 2.1: Define Nodes (Boxes, Text Labels, Fills, and Borders) ---------
define_flowchart_boxes <- function() {
  data.frame(
    xmin = c(
      25,                                # 1. Top source cohort
      2, 52,                             # 2-3. Empirical Arm | Simulation Arm
      2, 26,                             # 4-5. Strategy A | Strategy B
      2,                                 # 6. Prognostic models
      2,                                 # 7. Empirical outputs
      52,                                # 8. Reference population
      52,                                # 9. Missingness injection
      52, 61.4, 70.8, 80.2, 89.6,        # 10-14. S1, S2, S3, S3b, S4
      52                                 # 15. Simulation outputs
    ),
    xmax = c(
      75,
      48, 98,
      24, 48,
      48,
      48,
      98,
      98,
      60.4, 69.8, 79.2, 88.6, 98.0,
      98
    ),
    ymin = c(
      89,
      71, 71,
      54, 54,
      39,
      23,
      54,
      39,
      23, 23, 23, 23, 23,
      7
    ),
    ymax = c(
      97,
      84, 84,
      65, 65,
      49,
      33,
      65,
      49,
      33, 33, 33, 33, 33,
      17
    ),
    label = c(
      "Source Cohort: TCGA-LUAD & TCGA-LUSC\n(N = 867 primary NSCLC patients)",
      "Empirical Arm\n(Real-world data analysis)\nExpression (N = 867) | Methylation (N = 701)",
      "Simulation Arm\n(Monte Carlo benchmark with known truth)",
      "Strategy A:\nComplete-case\nexclusion",
      "Strategy B:\nUnknown-stage\nretention",
      "Prognostic models\n(LASSO, Elastic Net, Ridge, RSF, Stepwise Cox)",
      "Empirical outputs\nSelected genes, selection stability,\nC-index, methylation-score HR",
      "Reference population (N = 867)\nLatent severity for Unknown-stage patients;\nparameters estimated from the TCGA cohort",
      "Missingness injection (14.2%)\nMCAR | MAR | MNAR\n(survival times regenerated in each replicate)",
      "S1:\nComplete-case\ndeletion",
      "S2:\nUnknown\nretention",
      "S3:\nImputation\n(covariates)",
      "S3b:\nImputation\n(+ outcome)",
      "S4:\nIPW",
      "Simulation outputs\nBias, RMSE, coverage (100 independent replicates per mechanism)"
    ),
    fill = c(
      "#FFFFFF",
      "#FCF4F2", "#F0F5FA",
      "#FCF4F2", "#FCF4F2",
      "#FCF4F2",
      "#FCF4F2",
      "#F0F5FA",
      "#F0F5FA",
      "#F0F5FA", "#F0F5FA", "#F0F5FA", "#F0F5FA", "#F0F5FA",
      "#F0F5FA"
    ),
    color = rep("gray25", 15),
    size  = c(
      4.0,
      3.9, 3.9,
      3.6, 3.6,
      3.7,
      3.7,
      3.7,
      3.7,
      3.1, 3.1, 3.1, 3.1, 3.1,
      3.7
    ),
    stringsAsFactors = FALSE
  )
}

# --- Block 2.2: Define Directed Connectors (Arrows) ---------------------------
define_flowchart_arrows <- function() {
  data.frame(
    x    = c(50, 50, 25, 25, 13, 37, 25, 75, 75, 75,   75,   75,   75,   75,   56.2, 65.6, 75.0, 84.4, 93.8),
    y    = c(89, 89, 71, 71, 54, 54, 39, 71, 54, 39,   39,   39,   39,   39,   23,   23,   23,   23,   23),
    xend = c(25, 75, 13, 37, 25, 25, 25, 75, 75, 56.2, 65.6, 75.0, 84.4, 93.8, 75,   75,   75,   75,   75),
    yend = c(84, 84, 65, 65, 49, 49, 33, 65, 49, 33,   33,   33,   33,   33,   17,   17,   17,   17,   17)
  )
}


# ==============================================================================
# MODULE 3: Rendering & Export
# ==============================================================================

# --- Block 3.1: Build ggplot2 Flowchart & Save Publication Files (PDF/PNG) ----
build_study_flowchart <- function() {
  boxes  <- define_flowchart_boxes()
  arrows <- define_flowchart_arrows()
  
  ggplot() +
    geom_rect(
      data      = boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, color = color),
      linewidth = 0.75
    ) +
    geom_text(
      data       = boxes,
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label, size = size),
      family     = "serif",
      fontface   = "bold",
      lineheight = 1.12,
      color      = "black"
    ) +
    geom_segment(
      data      = arrows,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow     = arrow(length = unit(0.16, "cm"), type = "closed"),
      linewidth = 0.75,
      color     = "gray25"
    ) +
    scale_size_identity() +
    scale_fill_identity() +
    scale_color_identity() +
    coord_cartesian(xlim = c(0, 100), ylim = c(5, 99), expand = FALSE) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(10, 10, 10, 10)
    )
}

if (isTRUE(getOption("nsclc.fig01b_execute", TRUE))) {
  message("[Script 01b] Generating Figure 1 (Dual-arm study design flowchart)...")
  fig1 <- build_study_flowchart()
  
  ggsave("Figure1_Final.pdf", fig1, width = 12, height = 6.2, device = cairo_pdf)
  ggsave("Figure1_Final.png", fig1, width = 12, height = 6.2, dpi = 300, bg = "white")
  
  message("[Script 01b] Completed. Saved 'Figure1_Final.pdf' and 'Figure1_Final.png'.")
}
