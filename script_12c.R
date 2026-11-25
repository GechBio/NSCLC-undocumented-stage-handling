# ==============================================================================
# Script 12c: NSCLC Missing Data Empirical C-Index P-Values & Paired Comparisons
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Compute bootstrap standard errors, 95% confidence intervals, and
#           p-values for validation C-indices in the expression cohort:
#             (1) Discrimination vs. chance (C = 0.5) via Wald z-test
#             (2) Strategy B minus Strategy A for identical models (unpaired,
#                 descriptive across distinct validation sets)
#             (3) Gene/stepwise models vs. clinical-only baseline within each
#                 strategy via paired bootstrap on identical validation sets
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Base & Expression Model Objects
#   MODULE 2: Risk Prediction & Paired Bootstrap Engine
#       - Block 2.1: Harrell's C-Index & Model Risk Score Extractor
#       - Block 2.2: Simultaneous Paired Bootstrap Resampler (B = 2,000)
#   MODULE 3: Hypothesis Testing & Contrast Tables
#       - Block 3.1: Test (1) - Validation C-Index vs. 0.5 (Wald z)
#       - Block 3.2: Test (2) - Unpaired Strategy B minus Strategy A
#       - Block 3.3: Test (3) - Paired Model Contrast vs. Clinical Baseline
#   MODULE 4: Execution & Export
#       - Block 4.1: Run Full Pipeline & Save Supplementary CSV Tables (Table S3)
#
# Save as:  12c_cindex_hypothesis_tests.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_base.rds, nsclc_empirical_expr.rds
# Output:   TableS_Cindex_p_vs_0.5.csv, TableS_Cindex_B_minus_A.csv,
#           TableS_Cindex_vs_clinical_paired.csv
#
# Dependencies: glmnet, survival
#
# REVISION NOTES (Nov 2026):
#   - Refactored into standardized Module -> Block architecture matching Script 12.
#   - Added input file validation, explicit row.names = NULL formatting, and
#     execution guard (nsclc.emp12c_execute) for safe master-pipeline sourcing.
# ==============================================================================


# ==============================================================================
# MODULE 1: Environment & Dependency Management
# ==============================================================================

# --- Block 1.1: Load Required Packages ----------------------------------------
suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
})

# --- Block 1.2: Load Base & Expression Model Objects --------------------------
load_cindex_test_inputs <- function(
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
# MODULE 2: Risk Prediction & Paired Bootstrap Engine
# ==============================================================================

# --- Block 2.1: Harrell's C-Index & Model Risk Score Extractor ----------------
calc_cindex <- function(time, event, risk) {
  1 - survival::concordance(Surv(time, event) ~ risk)$concordance
}

extract_validation_risks <- function(base, emp, strategy) {
  genes   <- base$present_genes
  test_df <- base[[paste0("test_", strategy, "_expr")]]
  R       <- emp[[paste0("results_", strategy)]]
  x_mat   <- as.matrix(test_df[, genes])
  
  list(
    test = test_df,
    risk = list(
      clinical = as.numeric(predict(R$clinical$model, newdata = test_df, type = "lp")),
      lasso    = as.numeric(predict(R$lasso$model, newx = x_mat, s = "lambda.min")),
      enet     = as.numeric(predict(R$enet$model,  newx = x_mat, s = "lambda.min")),
      ridge    = as.numeric(predict(R$ridge$model, newx = x_mat, s = "lambda.min")),
      stepwise = as.numeric(predict(R$stepwise$model, newdata = test_df, type = "lp"))
    )
  )
}

# --- Block 2.2: Simultaneous Paired Bootstrap Resampler (B = 2,000) -----------
# Uses identical resample indices across all models within a strategy so
# within-strategy model differences are evaluated as paired bootstrap contrasts.
bootstrap_strategy_models <- function(base, emp, strategy, B = 2000, seed = 1) {
  g <- extract_validation_risks(base, emp, strategy)
  n <- nrow(g$test)
  set.seed(seed)
  
  est <- sapply(g$risk, function(r) {
    calc_cindex(g$test$os_time, g$test$os_event, r)
  })
  
  bs <- t(replicate(B, {
    idx <- sample.int(n, replace = TRUE)
    sapply(g$risk, function(r) {
      calc_cindex(g$test$os_time[idx], g$test$os_event[idx], r[idx])
    })
  }))
  
  list(est = est, bs = bs)
}


# ==============================================================================
# MODULE 3: Hypothesis Testing & Contrast Tables
# ==============================================================================

# --- Block 3.1: Test (1) - Validation C-Index vs. 0.5 (Wald z) ----------------
test_cindex_vs_chance <- function(boot_res) {
  do.call(rbind, lapply(c("A", "B"), function(S) {
    est <- boot_res[[S]]$est
    se  <- apply(boot_res[[S]]$bs, 2, sd, na.rm = TRUE)
    z   <- (est - 0.5) / se
    
    data.frame(
      Strategy = S,
      Model    = names(est),
      Cindex   = as.numeric(est),
      SE       = as.numeric(se),
      z        = as.numeric(z),
      p_wald   = 2 * pnorm(-abs(as.numeric(z))),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
}

# --- Block 3.2: Test (2) - Unpaired Strategy B minus Strategy A ---------------
test_strategy_b_minus_a <- function(boot_res) {
  mods <- names(boot_res$A$est)
  diff_est <- as.numeric(boot_res$B$est - boot_res$A$est)
  se_diff  <- as.numeric(sqrt(
    apply(boot_res$A$bs, 2, var, na.rm = TRUE) +
      apply(boot_res$B$bs, 2, var, na.rm = TRUE)
  ))
  z_val <- diff_est / se_diff
  
  data.frame(
    Model          = mods,
    Diff_B_minus_A = diff_est,
    SE             = se_diff,
    z              = z_val,
    p              = 2 * pnorm(-abs(z_val)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# --- Block 3.3: Test (3) - Paired Model Contrast vs. Clinical Baseline --------
test_paired_vs_clinical <- function(boot_res) {
  do.call(rbind, lapply(c("A", "B"), function(S) {
    bs  <- boot_res[[S]]$bs
    est <- boot_res[[S]]$est
    comp_models <- setdiff(colnames(bs), "clinical")
    
    do.call(rbind, lapply(comp_models, function(m) {
      d <- bs[, m] - bs[, "clinical"]
      data.frame(
        Strategy         = S,
        Model            = m,
        Diff_vs_clinical = as.numeric(est[m] - est["clinical"]),
        lo               = as.numeric(quantile(d, 0.025, na.rm = TRUE)),
        hi               = as.numeric(quantile(d, 0.975, na.rm = TRUE)),
        p_boot           = min(1, 2 * min(mean(d <= 0, na.rm = TRUE), mean(d >= 0, na.rm = TRUE))),
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }))
  }))
}


# ==============================================================================
# MODULE 4: Execution & Export
# ==============================================================================

# --- Block 4.1: Run Full Pipeline & Save Supplementary CSV Tables (Table S3) --
if (isTRUE(getOption("nsclc.emp12c_execute", TRUE))) {
  inputs <- load_cindex_test_inputs(
    base_path = "nsclc_base.rds",
    expr_path = "nsclc_empirical_expr.rds"
  )
  
  message("[Script 12c] Running paired bootstrap resamples (B = 2,000) across Strategies A & B...")
  boot_res <- list(
    A = bootstrap_strategy_models(inputs$base, inputs$emp, strategy = "A", B = 2000, seed = 1),
    B = bootstrap_strategy_models(inputs$base, inputs$emp, strategy = "B", B = 2000, seed = 1)
  )
  
  # (1) C-index vs 0.5
  t1 <- test_cindex_vs_chance(boot_res)
  cat("\n--- Test (1): Validation C-index vs. 0.5 (Wald z) ---\n")
  print(t1, digits = 3, row.names = FALSE)
  
  # (2) Strategy B minus Strategy A (unpaired descriptive)
  t2 <- test_strategy_b_minus_a(boot_res)
  cat("\n--- Test (2): Strategy B minus A (Unpaired, Descriptive) ---\n")
  print(t2, digits = 3, row.names = FALSE)
  
  # (3) Models vs. clinical-only baseline (paired bootstrap)
  t3 <- test_paired_vs_clinical(boot_res)
  cat("\n--- Test (3): Models vs. Clinical-Only Baseline (Paired Bootstrap) ---\n")
  print(t3, digits = 3, row.names = FALSE)
  
  write.csv(t1, file = "TableS_Cindex_p_vs_0.5.csv", row.names = FALSE)
  write.csv(t2, file = "TableS_Cindex_B_minus_A.csv", row.names = FALSE)
  write.csv(t3, file = "TableS_Cindex_vs_clinical_paired.csv", row.names = FALSE)
  
  message("[Script 12c] Completed. Saved TableS_Cindex_p_vs_0.5.csv, TableS_Cindex_B_minus_A.csv, and TableS_Cindex_vs_clinical_paired.csv.")
}
