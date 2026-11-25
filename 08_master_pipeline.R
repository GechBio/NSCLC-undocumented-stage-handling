# ==============================================================================
# Script 08: Master Execution Pipeline
#
# Author:   Getachew Derbew Belay
# Platform: R 4.5.2
#
# Purpose:  Run the full analysis sequentially and reproduce all manuscript
#           results, tables, figures and supplementary tables (S1-S5).
#
# ORDER MATTERS. Script 12 creates nsclc_repeated_splits.rds, which Script 03b
# needs, so the empirical sensitivity block (12-12e) runs BEFORE the figures.
# Helper files (05_simulation_monte_carlo.R and simulation_helpers.R) are not
# run directly as steps; they are sourced by 05b and 05c.
# ==============================================================================

steps <- c(
  # Step 1: base setup and data assembly (N = 867 expression, N = 701 methylation)
  "00_base_setup.R",
  # Step 2: empirical arm modelling
  "01_empirical_expression.R",
  "02_empirical_methylation.R",
  # Step 3: empirical sensitivity and inferential analyses (Tables S1-S4)
  "09_clinical_cohort_disambiguation.R",
  "script12_additional.R",              # Script 12  (bootstrap CIs, 100 splits, selection frequency)
  "12b_script_additional.R",            # Script 12b (bootstrap CIs: RSF, stepwise, methylation)
  "script_12c.R",                       # Script 12c (C-index p-values, paired comparisons)
  "Script_12d_Inferential_.R",          # Script 12d (Unknown category value, shared-split test)
  "12e_rsf_pvalues.R",                  # Script 12e (RSF p-values)
  # Step 4: empirical figures (Figures 2, 3 and 4)
  "03_generate_empirical_figures.R",
  "03b_feature_selection_upset.R",
  # Step 5: simulation arm (calibration, Table 4, Figures 5 and 6)
  "04_simulation_ground_truth.R",
  "04b_generate_simulation_ground_truth_figure.R",
  "05b_simulation_monte_carlo_independent.R",
  "06b_script_additional.R",            # Script 06b (Figure 6)
  # Step 6: stage-dependent MNAR scenarios (Table S5)
  "05c_simulation_stage_dependent_mnar.R",
  # Step 7: reporting checklist (Table 5)
  "07_reporting_checklist.R"
)
helpers <- c("simulation_helpers.R", "05_simulation_monte_carlo.R")

missing_files <- setdiff(c(steps, helpers), list.files())
if (length(missing_files) > 0L) {
  stop("These files were not found in the working directory:\n  ",
       paste(missing_files, collapse = "\n  "))
}

cat("========================================================\n")
cat("Starting full reproduction pipeline\n")
cat("========================================================\n")
for (f in steps) {
  cat("\n>>> EXECUTING:", f, "\n")
  source(f)
}
cat("========================================================\n")
cat("Pipeline complete: tables, figures and supplementary tables generated.\n")
cat("========================================================\n")
