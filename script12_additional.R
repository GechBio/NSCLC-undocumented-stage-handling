# ==============================================================================
# Script 12: NSCLC Missing Data Empirical Sensitivity & Stability Analyses
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Perform reporting-level empirical sensitivity and stability analyses
#           without altering primary modeling pipelines:
#             (a) Bootstrap 95% CIs of validation C-index (fixed models, B = 2,000)
#             (b) 100 repeated stratified 70/30 splits for C-index spread and
#                 LASSO gene selection frequency by strategy
#             (c) Unpenalized clinical-covariate forced-in LASSO Cox models
#             (d) Descriptive characterization & retention audit of Unknown stage
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Base & Empirical Model Objects
#   MODULE 2: Bootstrap Validation C-Index Confidence Intervals (Part a)
#       - Block 2.1: Concordance & Percentile Bootstrap Helpers
#       - Block 2.2: Evaluate Fixed Models Across Strategies A and B
#   MODULE 3: Repeated Stratified Splits & Selection Stability (Part b)
#       - Block 3.1: Single Stratified 70/30 Split Evaluator
#       - Block 3.2: 100-Split Execution & Frequency Summarization
#   MODULE 4: Forced-Clinical LASSO & Unknown-Stage Audit (Parts c & d)
#       - Block 4.1: Unpenalized Clinical Covariates LASSO Fit
#       - Block 4.2: Unknown-Stage Group Characterization & Retention Audit
#   MODULE 5: Execution & Export
#       - Block 5.1: Run Full Pipeline & Save Summary Tables / RDS Objects
#
# Save as:  12_empirical_sensitivity_analyses.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_base.rds, nsclc_empirical_expr.rds
# Output:   TableS_Cindex_bootstrap_CI.csv, nsclc_repeated_splits.rds
#
# Dependencies: glmnet, survival, caret
#
# REVISION NOTES (Nov 2026):
#   - Refactored into standardized Module -> Block architecture matching Script 05.
#   - Preserved numeric column types in TableS_Cindex_bootstrap_CI.csv by
#     constructing rows with data.frame() rather than character-coercing c().
#   - Added input file validation and execution guard (nsclc.emp12_execute).
# ==============================================================================


# ==============================================================================
# MODULE 1: Environment & Dependency Management
# ==============================================================================

# --- Block 1.1: Load Required Packages ----------------------------------------
suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(caret)
})

# --- Block 1.2: Load Base & Empirical Model Objects ---------------------------
load_empirical_inputs <- function(
    base_path = "nsclc_base.rds",
    emp_path  = "nsclc_empirical_expr.rds"
) {
  if (!file.exists(base_path)) stop("Required input file not found: ", base_path)
  if (!file.exists(emp_path))  stop("Required input file not found: ", emp_path)
  
  list(
    base = readRDS(base_path),
    emp  = readRDS(emp_path)
  )
}


# ==============================================================================
# MODULE 2: Bootstrap Validation C-Index Confidence Intervals (Part a)
# ==============================================================================

# --- Block 2.1: Concordance & Percentile Bootstrap Helpers --------------------
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

# --- Block 2.2: Evaluate Fixed Models Across Strategies A and B ---------------
run_bootstrap_cindex_table <- function(base, emp, B = 2000, seed = 1) {
  genes <- base$present_genes
  res   <- list()
  
  for (S in c("A", "B")) {
    test_df <- base[[paste0("test_", S, "_expr")]]
    mod_res <- emp[[paste0("results_", S)]]
    
    for (m in c("lasso", "enet", "ridge")) {
      lp <- as.numeric(predict(
        mod_res[[m]]$model,
        newx = as.matrix(test_df[, genes]),
        s    = "lambda.min"
      ))
      row_df <- cbind(
        data.frame(Strategy = S, Model = m, stringsAsFactors = FALSE),
        boot_cindex_ci(test_df, lp, B = B, seed = seed)
      )
      res[[paste(S, m, sep = "_")]] <- row_df
    }
    
    lp_clin <- as.numeric(predict(mod_res$clinical$model, newdata = test_df, type = "lp"))
    row_clin <- cbind(
      data.frame(Strategy = S, Model = "clinical", stringsAsFactors = FALSE),
      boot_cindex_ci(test_df, lp_clin, B = B, seed = seed)
    )
    res[[paste(S, "clinical", sep = "_")]] <- row_clin
  }
  
  do.call(rbind, res)
}


# ==============================================================================
# MODULE 3: Repeated Stratified Splits & Selection Stability (Part b)
# ==============================================================================

# --- Block 3.1: Single Stratified 70/30 Split Evaluator -------------------------
evaluate_one_split <- function(cohort, genes, seed) {
  set.seed(seed)
  strata <- interaction(cohort$os_event, cohort$stage_group, cohort$histology, drop = TRUE)
  idx    <- caret::createDataPartition(strata, p = 0.7, list = FALSE)
  
  tr <- cohort[idx, ]
  te <- cohort[-idx, ]
  
  cv_fit <- cv.glmnet(
    x            = as.matrix(tr[, genes]),
    y            = Surv(tr$os_time, tr$os_event),
    family       = "cox",
    alpha        = 1,
    nfolds       = 10,
    type.measure = "C"
  )
  
  lp_lasso <- as.numeric(predict(cv_fit, newx = as.matrix(te[, genes]), s = "lambda.min"))
  cf       <- coef(cv_fit, s = "lambda.min")
  sel      <- rownames(cf)[as.numeric(cf) != 0]
  
  clin_fit <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = tr)
  lp_clin  <- as.numeric(predict(clin_fit, newdata = te, type = "lp"))
  
  list(
    c      = calc_cindex(te$os_time, te$os_event, lp_lasso),
    c_clin = calc_cindex(te$os_time, te$os_event, lp_clin),
    sel    = sel
  )
}

# --- Block 3.2: 100-Split Execution & Frequency Summarization -----------------
run_repeated_splits <- function(base, n_splits = 100, base_seed = 1000) {
  genes   <- base$present_genes
  rep_out <- list()
  
  for (S in c("A", "B")) {
    coh <- base[[paste0("cohort_", S, "_expr")]]
    r   <- lapply(seq_len(n_splits), function(k) {
      evaluate_one_split(coh, genes = genes, seed = base_seed + k)
    })
    
    c_vec      <- sapply(r, `[[`, "c")
    c_clin_vec <- sapply(r, `[[`, "c_clin")
    n_sel_vec  <- sapply(r, function(z) length(z$sel))
    freq_tab   <- sort(
      table(factor(unlist(lapply(r, `[[`, "sel")), levels = genes)) / n_splits,
      decreasing = TRUE
    )
    
    rep_out[[S]] <- list(
      c      = c_vec,
      c_clin = c_clin_vec,
      n_sel  = n_sel_vec,
      freq   = freq_tab
    )
    
    cat(
      "\nStrategy", S,
      ": LASSO C-index median", round(median(c_vec), 3),
      "[2.5-97.5%]", round(quantile(c_vec, c(0.025, 0.975)), 3),
      "| clinical median", round(median(c_clin_vec), 3),
      "| median genes selected", median(n_sel_vec), "\n"
    )
    print(round(freq_tab, 2))
  }
  
  return(rep_out)
}


# ==============================================================================
# MODULE 4: Forced-Clinical LASSO & Unknown-Stage Audit (Parts c & d)
# ==============================================================================

# --- Block 4.1: Unpenalized Clinical Covariates LASSO Fit ---------------------
fit_forced_clinical_lasso <- function(train_df, genes, seed = 123) {
  mm <- model.matrix(~ age + gender + histology + stage_group, data = train_df)[, -1, drop = FALSE]
  x  <- cbind(mm, as.matrix(train_df[, genes]))
  pf <- c(rep(0, ncol(mm)), rep(1, length(genes)))
  
  set.seed(seed)
  cv_fit <- cv.glmnet(
    x              = x,
    y              = Surv(train_df$os_time, train_df$os_event),
    family         = "cox",
    alpha          = 1,
    nfolds         = 10,
    type.measure   = "C",
    penalty.factor = pf
  )
  
  cf <- coef(cv_fit, s = "lambda.min")
  setdiff(rownames(cf)[as.numeric(cf) != 0], colnames(mm))
}

# --- Block 4.2: Unknown-Stage Group Characterization & Retention Audit --------
audit_unknown_stage <- function(base) {
  B   <- base$cohort_B_expr
  cov <- base$covariates_clean
  
  cat("\n--- Part (d): Unknown-Stage Group Characterization ---\n")
  print(aggregate(cbind(os_event, os_time, age) ~ stage_group, data = B, FUN = mean))
  print(prop.table(table(B$stage_group, B$histology), margin = 1))
  print(table(B$stage_group, substr(B$patient_id, 6, 7)))  # Tissue source site (TSS)
  
  cat(
    "\nUnknown share before survival/NA filters:", round(mean(cov$stage_group == "Unknown"), 4),
    "| after filters:", round(mean(B$stage_group == "Unknown"), 4), "\n"
  )
  cat("Retained fraction by stage (clinical covariate file -> final analytic cohort):\n")
  print(round(tapply(cov$patient_id %in% B$patient_id, cov$stage_group, mean), 4))
}


# ==============================================================================
# MODULE 5: Execution & Export
# ==============================================================================

# --- Block 5.1: Run Full Pipeline & Save Summary Tables / RDS Objects ---------
if (isTRUE(getOption("nsclc.emp12_execute", TRUE))) {
  inputs <- load_empirical_inputs("nsclc_base.rds", "nsclc_empirical_expr.rds")
  
  base <- inputs$base
  emp  <- inputs$emp

  cat("\n--- Part (a): Bootstrap 95% CIs (B = 2,000) ---\n")
  boot_tab <- run_bootstrap_cindex_table(base, emp, B = 2000, seed = 1)
  print(boot_tab)
  write.csv(boot_tab, "TableS_Cindex_bootstrap_CI.csv", row.names = FALSE)

  cat("\n--- Part (b): 100 Repeated Stratified Splits ---\n")
  rep_splits <- run_repeated_splits(base, n_splits = 100, base_seed = 1000)
  saveRDS(rep_splits, "nsclc_repeated_splits.rds")

  cat("\n--- Part (c): Forced-Clinical LASSO Selected Genes ---\n")
  genes <- base$present_genes
  sel_forced_A <- fit_forced_clinical_lasso(base$train_A_expr, genes, seed = 123)
  sel_forced_B <- fit_forced_clinical_lasso(base$train_B_expr, genes, seed = 123)
  cat("Strategy A forced-clinical LASSO selected:", paste(sel_forced_A, collapse = ", "), "\n")
  cat("Strategy B forced-clinical LASSO selected:", paste(sel_forced_B, collapse = ", "), "\n")

  audit_unknown_stage(base)

  cat("\nSaved TableS_Cindex_bootstrap_CI.csv and nsclc_repeated_splits.rds.\n")
  cat("Script 12 Complete!\n")
}
