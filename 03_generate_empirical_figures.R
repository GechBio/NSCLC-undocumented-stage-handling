# ==============================================================================
# Script 03: NSCLC Missing Data Visualization - Empirical Arm Figures (Figures 2 & 4)
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Generate empirical arm publication figures in R (ggplot2):
#             - Figure 2: Validation C-index with 95% percentile bootstrap CIs
#               and right-hand CI table across six models (Strategy A vs. B)
#             - Figure 4: Annotated forest plot of adjusted hazard ratio per IQR
#               for the gene-level methylation risk score
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Empirical Results & Bootstrap CI Tables
#   MODULE 2: Figure 2 - Validation C-Index Comparison (Expression Cohort)
#       - Block 2.1: Assemble Point Estimates & 95% Bootstrap CIs
#       - Block 2.2: Render & Export Figure 2 (PDF & 300-DPI PNG)
#   MODULE 3: Figure 4 - Annotated Forest Plot (Methylation Cohort)
#       - Block 3.1: Format Adjusted Hazard Ratios & Non-Overlapping Labels
#       - Block 3.2: Render & Export Figure 4 (PDF & 300-DPI PNG)
#   MODULE 4: Execution & Export
#       - Block 4.1: Run Figure Generation Pipeline
#
# Save as:  03_generate_empirical_figures.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_empirical_expr.rds, nsclc_empirical_meth.rds,
#           TableS_Cindex_bootstrap_CI.csv, TableS_Cindex_bootstrap_CI_extra.csv
# Output:   Figure2_Cindex_Comparison.pdf (.png)
#           Figure4_Methylation_HR_Annotated.pdf (.png)
#
# Dependencies: ggplot2, dplyr, tidyr
# ==============================================================================


# ==============================================================================
# MODULE 1: Environment & Dependency Management
# ==============================================================================

# --- Block 1.1: Load Required Packages ----------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# --- Block 1.2: Load Empirical Results & Bootstrap CI Tables ------------------
load_empirical_figure_inputs <- function(
    expr_path    = "nsclc_empirical_expr.rds",
    meth_path    = "nsclc_empirical_meth.rds",
    boot12_path  = "TableS_Cindex_bootstrap_CI.csv",
    boot12b_path = "TableS_Cindex_bootstrap_CI_extra.csv"
) {
  if (!file.exists(expr_path) || !file.exists(meth_path)) {
    stop("Run Scripts 01 and 02 first to generate nsclc_empirical_expr.rds and nsclc_empirical_meth.rds.")
  }
  
  has_boot <- file.exists(boot12_path) && file.exists(boot12b_path)
  if (!has_boot) {
    warning("Bootstrap CSVs from Scripts 12/12b not found; Figure 2 will plot point estimates only.")
  }
  
  list(
    expr     = readRDS(expr_path),
    meth     = readRDS(meth_path),
    has_boot = has_boot,
    boot12   = if (has_boot) read.csv(boot12_path,  stringsAsFactors = FALSE) else NULL,
    boot12b  = if (has_boot) read.csv(boot12b_path, stringsAsFactors = FALSE) else NULL
  )
}


# ==============================================================================
# MODULE 2: Figure 2 - Validation C-Index Comparison (Expression Cohort)
# ==============================================================================

# --- Block 2.1: Assemble Point Estimates & 95% Bootstrap CIs ------------------
prepare_fig2_data <- function(inputs) {
  model_levels <- c("Ridge", "LASSO", "Elastic Net", "RSF", "Stepwise Cox", "Clinical only")
  
  if (inputs$has_boot) {
    b1 <- inputs$boot12 %>%
      filter(Model %in% c("clinical", "lasso", "enet", "ridge")) %>%
      mutate(Model_Label = recode(
        Model,
        "clinical" = "Clinical only",
        "lasso"    = "LASSO",
        "enet"     = "Elastic Net",
        "ridge"    = "Ridge"
      ))
    
    b2 <- inputs$boot12b %>%
      filter(Model %in% c("expr stepwise", "expr rsf")) %>%
      mutate(Model_Label = recode(
        Model,
        "expr stepwise" = "Stepwise Cox",
        "expr rsf"      = "RSF"
      ))
    
    df <- bind_rows(b1, b2) %>% dplyr::rename(Cindex = 3, CI_lower = 4, CI_upper = 5) %>%
      mutate(
        Model_Label = factor(Model_Label, levels = model_levels),
        Strategy_Label = ifelse(
          Strategy == "A",
          "Strategy A (complete-case, N = 744)",
          "Strategy B (Unknown retained, N = 867)"
        ),
        CI_Text = sprintf("%.3f (%.3f–%.3f)", Cindex, CI_lower, CI_upper),
        Y_Num   = as.numeric(Model_Label) + ifelse(Strategy == "B", 0.14, -0.14),
        has_ci  = TRUE
      )
    return(df)
  }
  
  # Fallback if bootstrap CSVs have not been run yet
  expr_data <- inputs$expr
  cindex_df <- data.frame(
    Model_Label = c("Clinical only", "Stepwise Cox", "RSF", "Elastic Net", "LASSO", "Ridge"),
    Strategy_A  = expr_data$cindex_table$Strategy_A_Cindex,
    Strategy_B  = expr_data$cindex_table$Strategy_B_Cindex,
    stringsAsFactors = FALSE
  )
  
  cindex_df %>%
    pivot_longer(
      cols      = c("Strategy_A", "Strategy_B"),
      names_to  = "Strategy",
      values_to = "Cindex"
    ) %>%
    mutate(
      Model_Label = factor(Model_Label, levels = model_levels),
      Strategy_Label = ifelse(
        Strategy == "Strategy_A",
        "Strategy A (complete-case, N = 744)",
        "Strategy B (Unknown retained, N = 867)"
      ),
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      CI_Text  = sprintf("%.3f", Cindex),
      Y_Num    = as.numeric(Model_Label) + ifelse(Strategy == "Strategy_B", 0.14, -0.14),
      has_ci   = FALSE
    )
}

# --- Block 2.2: Render & Export Figure 2 (PDF & 300-DPI PNG) ------------------
build_figure2 <- function(fig2_df) {
  pos <- position_dodge(width = 0.42)
  
  p <- ggplot(fig2_df, aes(x = Cindex, y = Model_Label, color = Strategy_Label)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray55", linewidth = 0.7)
  
  if (all(fig2_df$has_ci)) {
    p <- p +
      geom_errorbarh(
        aes(xmin = CI_lower, xmax = CI_upper),
        height    = 0,
        linewidth = 1.1,
        position  = pos
      ) +
      geom_point(size = 3.5, position = pos)
  } else {
    p <- p +
      geom_line(aes(group = Model_Label), color = "gray60", linewidth = 1) +
      geom_point(size = 3.5, position = pos)
  }
  
  p +
    annotate(
      "text",
      x = 0.755, y = 6.55,
      label = "C-index (95% CI)",
      hjust = 0, fontface = "bold", size = 3.8, color = "black"
    ) +
    geom_text(
      aes(x = 0.755, y = Y_Num, label = CI_Text, color = Strategy_Label),
      hjust = 0, size = 3.5, show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      "Strategy A (complete-case, N = 744)"    = "#D55E00",
      "Strategy B (Unknown retained, N = 867)" = "#0072B2"
    )) +
    scale_x_continuous(
      limits = c(0.38, 0.89),
      breaks = seq(0.40, 0.70, by = 0.05)
    ) +
    coord_cartesian(ylim = c(0.6, 6.6), clip = "off") +
    labs(
      title   = "Validation C-index with 95% confidence intervals (expression cohort)",
      caption = "Percentile bootstrap CIs (2,000 resamples of the validation set, models fixed).",
      x       = "Validation C-index (95% CI)",
      y       = NULL,
      color   = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.text     = element_text(size = 10.5),
      plot.title      = element_text(face = "bold", size = 13, hjust = 0),
      plot.caption    = element_text(size = 9, hjust = 0, color = "gray30", margin = margin(t = 8)),
      axis.text.y     = element_text(face = "bold", size = 11, color = "black"),
      axis.text.x     = element_text(size = 10, color = "black"),
      axis.title.x    = element_text(size = 11, margin = margin(t = 8)),
      plot.margin     = margin(12, 20, 10, 10)
    )
}


# ==============================================================================
# MODULE 3: Figure 4 - Annotated Forest Plot (Methylation Cohort)
# ==============================================================================

# --- Block 3.1: Format Adjusted Hazard Ratios & Non-Overlapping Labels --------
prepare_fig4_data <- function(meth_data) {
  hr_df <- data.frame(
    Strategy = factor(
      c("Strategy A\n(Complete-Case)", "Strategy B\n(Missing-Retained)"),
      levels = c("Strategy A\n(Complete-Case)", "Strategy B\n(Missing-Retained)")
    ),
    HR    = c(meth_data$results_A$gene_level$hr_iqr,    meth_data$results_B$gene_level$hr_iqr),
    Lower = c(meth_data$results_A$gene_level$hr_ci[1], meth_data$results_B$gene_level$hr_ci[1]),
    Upper = c(meth_data$results_A$gene_level$hr_ci[2], meth_data$results_B$gene_level$hr_ci[2]),
    Pval  = c(meth_data$results_A$gene_level$p_value,   meth_data$results_B$gene_level$p_value)
  )
  
  hr_df$Label <- sprintf(
    "HR: %.2f (95%% CI: %.2f–%.2f)\np = %.3f",
    hr_df$HR, hr_df$Lower, hr_df$Upper, hr_df$Pval
  )
  
  # Shift Strategy B label strictly right of x = 1.00 to avoid colliding with reference line
  hr_df$Label_X     <- ifelse(grepl("Strategy B", hr_df$Strategy), 1.03, hr_df$HR)
  hr_df$Label_Hjust <- ifelse(grepl("Strategy B", hr_df$Strategy), 0,    0.5)
  
  return(hr_df)
}

# --- Block 3.2: Render & Export Figure 4 (PDF & 300-DPI PNG) ------------------
build_figure4 <- function(hr_df) {
  ggplot(hr_df, aes(x = HR, y = Strategy, color = Strategy)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_pointrange(aes(xmin = Lower, xmax = Upper), size = 1.1, linewidth = 1.1) +
    geom_text(
      aes(x = Label_X, label = Label, hjust = Label_Hjust),
      vjust    = -0.95,
      size     = 4.0,
      fontface = "bold",
      color    = "black"
    ) +
    scale_color_manual(values = c(
      "Strategy A\n(Complete-Case)"    = "#D55E00",
      "Strategy B\n(Missing-Retained)" = "#0072B2"
    )) +
    scale_x_continuous(limits = c(0.8, 1.9), breaks = seq(0.8, 1.8, by = 0.2)) +
    scale_y_discrete(expand = expansion(mult = c(0.35, 0.55))) +
    labs(
      title    = "Methylation Risk Score: Adjusted HR per IQR",
      subtitle = "Adjusted for Age, Gender, Histology, and Stage (Strategy B: N = 701; Strategy A: N = 603)",
      x        = "Hazard Ratio (95% Confidence Interval)",
      y        = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title      = element_text(face = "bold", size = 13.5, hjust = 0.5),
      plot.subtitle   = element_text(size = 10.5, hjust = 0.5, margin = margin(b = 12)),
      axis.text.y     = element_text(face = "bold", size = 11, color = "black"),
      axis.text.x     = element_text(size = 10.5, color = "black"),
      axis.title.x    = element_text(face = "bold", size = 11.5, margin = margin(t = 10))
    )
}


# ==============================================================================
# MODULE 4: Execution & Export
# ==============================================================================

# --- Block 4.1: Run Figure Generation Pipeline --------------------------------
if (isTRUE(getOption("nsclc.fig03_execute", TRUE))) {
  inputs <- load_empirical_figure_inputs()
  
  message("[Script 03] Generating Figure 2 (Validation C-index comparison)...")
  fig2_df <- prepare_fig2_data(inputs)
  fig2    <- build_figure2(fig2_df)
  ggsave("Figure2_Cindex_Comparison.pdf", fig2, width = 9.2, height = 5.2, device = cairo_pdf)
  ggsave("Figure2_Cindex_Comparison.png", fig2, width = 9.2, height = 5.2, dpi = 300, bg = "white")
  
  message("[Script 03] Generating Figure 4 (Methylation HR forest plot)...")
  fig4_df <- prepare_fig4_data(inputs$meth)
  fig4    <- build_figure4(fig4_df)
  ggsave("Figure4_Methylation_HR_Annotated.pdf", fig4, width = 8.0, height = 5.0, device = cairo_pdf)
  ggsave("Figure4_Methylation_HR_Annotated.png", fig4, width = 8.0, height = 5.0, dpi = 300, bg = "white")
  
  message("[Script 03] Completed. Saved Figure 2 and Figure 4 (PDF and 300-DPI PNG).")
}
