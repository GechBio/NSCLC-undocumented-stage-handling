# ==============================================================================
# Script 12e: NSCLC Missing Data Exact Bootstrap P-Values for RSF Validation C-Index
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Compute exact bootstrap standard errors, Wald z-statistics, and
#           two-sided p-values (Wald z and percentile bootstrap) for Random
#           Survival Forest (RSF) validation C-index vs. chance (C = 0.5)
#           across Strategies A and B, completing Table 1 of the manuscript.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Base & Expression Model Objects
#   MODULE 2: RSF Discrimination & Bootstrap Evaluation
#       - Block 2.1: Harrell's C-Index & Bootstrap Hypothesis Test Helper
#   MODULE 3: Execution & Export
#       - Block 3.1: Evaluate Strategies A & B and Save Supplementary CSV Table
#
# Save as:  12e_rsf_cindex_pvalues.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_base.rds, nsclc_empirical_expr.rds
# Output:   TableS_RSF_p_vs_0.5.csv
#
# Dependencies: survival, randomForestSRC
#
# REVISION NOTES (Nov 2026):
#   - Refactored into standardized Module -> Block architecture matching Script 12.
#   - Added input file validation, explicit row.names = NULL formatting, and
#     execution guard (nsclc.emp12e_execute) for safe master-pipeline sourcing.
# ==============================================================================


# ==============================================================================
# MODULE 1: Environment & Dependency Management
# ==============================================================================

# --- Block 1.1: Load Required Packages ----------------------------------------
suppressPackageStartupMessages({
  library(survival)
  library(randomForestSRC)
})

# --- Block 1.2: Load Base & Expression Model Objects --------------------------
load_rsf_test_inputs <- function(
    base_path = "nsclc_base.rds",
    expr_path = "nsclc_empirical_expr.rds"
) {
  if (!file.exists(base_path)) stop("Required input file not found: ", base_path)
  if (!file.exists(expr_path)) stop("Required input file not found: ", expr_path)
  
  list(
    base = readRDS(base_path),
    emp  = readRDS(expr_path)
  )
}


# ==============================================================================
# MODULE 2: RSF Discrimination & Bootstrap Evaluation
# ==============================================================================

# --- Block 2.1: Harrell's C-Index & Bootstrap Hypothesis Test Helper ----------
calc_cindex <- function(time, event, risk) {
  1 - survival::concordance(Surv(time, event) ~ risk)$concordance
}

evaluate_rsf_vs_chance <- function(base, emp, strategy, B = 2000, seed = 1) {
  genes   <- base$present_genes
  test_df <- base[[paste0("test_", strategy, "_expr")]]
  R       <- emp[[paste0("results_", strategy)]]
  
  risk <- as.numeric(
    predict(R$rsf$model, newdata = test_df[, c("os_time", "os_event", genes)])$predicted
  )
  
  set.seed(seed)
  n <- nrow(test_df)
  boot_vals <- replicate(B, {
    idx <- sample.int(n, replace = TRUE)
    calc_cindex(test_df$os_time[idx], test_df$os_event[idx], risk[idx])
  })
  
  est    <- calc_cindex(test_df$os_time, test_df$os_event, risk)
  se     <- sd(boot_vals, na.rm = TRUE)
  z_val  <- (est - 0.5) / se
  p_wald <- 2 * pnorm(-abs(z_val))
  p_boot <- min(1, 2 * min(mean(boot_vals <= 0.5, na.rm = TRUE), mean(boot_vals >= 0.5, na.rm = TRUE)))
  
  data.frame(
    Strategy  = strategy,
    Model     = "rsf",
    Cindex    = est,
    SE        = se,
    z         = z_val,
    p_wald    = p_wald,
    p_boot    = p_boot,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# MODULE 3: Execution & Export
# ==============================================================================

# --- Block 3.1: Evaluate Strategies A & B and Save Supplementary CSV Table ----
run_rsf_pvalues <- function(base, emp, B = 2000, seed = 1) {
  do.call(rbind, lapply(c("A", "B"), function(S) {
    evaluate_rsf_vs_chance(base = base, emp = emp, strategy = S, B = B, seed = seed)
  }))
}

if (isTRUE(getOption("nsclc.emp12e_execute", TRUE))) {
  inputs <- load_rsf_test_inputs(
    base_path = "nsclc_base.rds",
    expr_path = "nsclc_empirical_expr.rds"
  )
  
  message("[Script 12e] Computing exact bootstrap p-values for RSF C-index vs. 0.5 (B = 2,000)...")
  tab_rsf <- run_rsf_pvalues(base = inputs$base, emp = inputs$emp, B = 2000, seed = 1)
  
  print(tab_rsf, digits = 3, row.names = FALSE)
  write.csv(tab_rsf, file = "TableS_RSF_p_vs_0.5.csv", row.names = FALSE)
  
  message("[Script 12e] Completed. Saved 'TableS_RSF_p_vs_0.5.csv'.")
}
