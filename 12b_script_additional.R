# ==============================================================================
# Script 12b: NSCLC Missing Data Empirical Bootstrap CIs - RSF, Stepwise & Methylation
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Compute 95% percentile bootstrap confidence intervals (B = 2,000)
#           for validation C-indices of Random Survival Forest (RSF), Stepwise
#           Cox, and methylation models (gene-level signature and clinical
#           baseline), completing Table 1 of the manuscript.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Base, Expression, and Methylation Model Objects
#   MODULE 2: Concordance & Bootstrap Helpers
#       - Block 2.1: Harrell's C-Index & Percentile Bootstrap Functions
#   MODULE 3: Model Evaluation Across Modalities (Strategies A & B)
#       - Block 3.1: Expression Models (Stepwise Cox & RSF)
#       - Block 3.2: Methylation Models (LASSO Signature & Clinical Baseline)
#   MODULE 4: Execution & Export
#       - Block 4.1: Combine Results & Save Supplementary CSV Table
#
# Save as:  12b_bootstrap_cindex_extra.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_base.rds, nsclc_empirical_expr.rds, nsclc_empirical_meth.rds
# Output:   TableS_Cindex_bootstrap_CI_extra.csv
#
# Dependencies: survival, glmnet, randomForestSRC
#
# REVISION NOTES (Nov 2026):
#   - Refactored into standardized Module -> Block architecture matching Script 12.
#   - Preserved numeric column types in TableS_Cindex_bootstrap_CI_extra.csv by
#     constructing rows with data.frame() instead of character-coercing c().
#   - Added zero-feature guard for methylation LASSO signature scoring and an
#     execution guard (nsclc.emp12b_execute) for safe pipeline sourcing.
# ==============================================================================


# ==============================================================================
# MODULE 1: Environment & Dependency Management
# ==============================================================================

# --- Block 1.1: Load Required Packages ----------------------------------------
suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(randomForestSRC)
})

# --- Block 1.2: Load Base, Expression, and Methylation Model Objects ----------
load_extra_bootstrap_inputs <- function(
    base_path = "nsclc_base.rds",
    expr_path = "nsclc_empirical_expr.rds",
    meth_path = "nsclc_empirical_meth.rds"
) {
  if (!file.exists(base_path)) stop("Required input file not found: ", base_path)
  if (!file.exists(expr_path)) stop("Required input file not found: ", expr_path)
  if (!file.exists(meth_path)) stop("Required input file not found: ", meth_path)
  
  list(
    base = readRDS(base_path),
    emp  = readRDS(expr_path),
    meth = readRDS(meth_path)
  )
}


# ==============================================================================
# MODULE 2: Concordance & Bootstrap Helpers
# ==============================================================================

# --- Block 2.1: Harrell's C-Index & Percentile Bootstrap Functions ------------
calc_cindex <- function(time, event, risk) {
  1 - survival::concordance(Surv(time, event) ~ risk)$concordance
}

boot_cindex_ci <- function(test_df, risk, B = 2000, seed = 1) {
  set.seed(seed)
  n <- nrow(test_df)
  
  boot_vals <- replicate(B, {
    idx <- sample.int(n, replace = TRUE)
    calc_cindex(test_df$os_time[idx], test_df$os_event[idx], risk[idx])
  })
  
  qs <- quantile(boot_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  data.frame(
    Cindex   = calc_cindex(test_df$os_time, test_df$os_event, risk),
    CI_lower = as.numeric(qs[1]),
    CI_upper = as.numeric(qs[2]),
    Boot_SE  = sd(boot_vals, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# MODULE 3: Model Evaluation Across Modalities (Strategies A & B)
# ==============================================================================

# --- Block 3.1: Expression Models (Stepwise Cox & RSF) ------------------------
evaluate_expr_extra <- function(base, emp, strategy, B = 2000, seed = 1) {
  genes   <- base$present_genes
  test_df <- base[[paste0("test_", strategy, "_expr")]]
  R       <- emp[[paste0("results_", strategy)]]
  
  # Stepwise Cox linear predictor
  lp_step <- as.numeric(predict(R$stepwise$model, newdata = test_df, type = "lp"))
  row_step <- cbind(
    data.frame(Strategy = strategy, Model = "expr stepwise", stringsAsFactors = FALSE),
    boot_cindex_ci(test_df, lp_step, B = B, seed = seed)
  )
  
  # Random Survival Forest predicted mortality
  pr_rsf <- predict(R$rsf$model, newdata = test_df[, c("os_time", "os_event", genes)])$predicted
  row_rsf <- cbind(
    data.frame(Strategy = strategy, Model = "expr rsf", stringsAsFactors = FALSE),
    boot_cindex_ci(test_df, as.numeric(pr_rsf), B = B, seed = seed)
  )
  
  rbind(row_step, row_rsf)
}

# --- Block 3.2: Methylation Models (LASSO Signature & Clinical Baseline) ------
evaluate_meth_extra <- function(base, meth, strategy, B = 2000, seed = 1) {
  trm <- base[[paste0("train_", strategy, "_meth")]]
  tem <- base[[paste0("test_", strategy, "_meth")]]
  M   <- meth[[paste0("results_", strategy)]]$gene_level
  
  # Gene-level LASSO methylation risk score
  cf   <- coef(M$model, s = "lambda.min")
  nz   <- which(as.numeric(cf) != 0)
  act  <- rownames(cf)[nz]
  beta <- as.numeric(cf)[nz]
  
  if (length(act) > 0L) {
    rs_meth <- as.numeric(as.matrix(tem[, act, drop = FALSE]) %*% beta)
  } else {
    rs_meth <- rep(0, nrow(tem))
  }
  
  row_sig <- cbind(
    data.frame(Strategy = strategy, Model = "meth signature", stringsAsFactors = FALSE),
    boot_cindex_ci(tem, rs_meth, B = B, seed = seed)
  )
  
  # Methylation-cohort clinical baseline Cox model
  clin_fit <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = trm)
  lp_clin  <- as.numeric(predict(clin_fit, newdata = tem, type = "lp"))
  
  row_clin <- cbind(
    data.frame(Strategy = strategy, Model = "meth clinical", stringsAsFactors = FALSE),
    boot_cindex_ci(tem, lp_clin, B = B, seed = seed)
  )
  
  rbind(row_sig, row_clin)
}


# ==============================================================================
# MODULE 4: Execution & Export
# ==============================================================================

# --- Block 4.1: Combine Results & Save Supplementary CSV Table ----------------
run_bootstrap_extra <- function(base, emp, meth, B = 2000, seed = 1) {
  res_list <- lapply(c("A", "B"), function(S) {
    rbind(
      evaluate_expr_extra(base, emp, strategy = S, B = B, seed = seed),
      evaluate_meth_extra(base, meth, strategy = S, B = B, seed = seed)
    )
  })
  do.call(rbind, res_list)
}

if (isTRUE(getOption("nsclc.emp12b_execute", TRUE))) {
  inputs <- load_extra_bootstrap_inputs(
    base_path = "nsclc_base.rds",
    expr_path = "nsclc_empirical_expr.rds",
    meth_path = "nsclc_empirical_meth.rds"
  )
  
  message("[Script 12b] Computing bootstrap 95% CIs for RSF, Stepwise Cox, and Methylation...")
  tab_extra <- run_bootstrap_extra(
    base = inputs$base,
    emp  = inputs$emp,
    meth = inputs$meth,
    B    = 2000,
    seed = 1
  )
  
  print(tab_extra, digits = 3, row.names = FALSE)
  write.csv(tab_extra, file = "TableS_Cindex_bootstrap_CI_extra.csv", row.names = FALSE)
  
  message("[Script 12b] Completed. Saved 'TableS_Cindex_bootstrap_CI_extra.csv'.")
}
