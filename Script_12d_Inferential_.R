# ==============================================================================
# Script 12d: NSCLC Missing Data Inferential Comparison of Strategy A vs. Strategy B
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Perform inferential tests comparing Strategy A (Complete-Case) vs.
#           Strategy B (Unknown-Stage Retained) and verify simulation inputs:
#             (0) Full-cohort calibration Cox model: stage HRs with vs. without
#                 full-cohort LASSO genes (verifies HR 3.22 and HR 1.86)
#             (1) Predictive value of the Unknown category on a shared validation
#                 set (Strategy B primary test split) via paired bootstrap delta
#                 C-index and full-cohort likelihood-ratio test
#             (2) Training-set utility of retaining Unknown stage across 100
#                 shared splits scored on identical documented-stage test sets,
#                 evaluated via Nadeau-Bengio corrected paired t-tests
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Base Cohort & Split Objects
#   MODULE 2: Full-Cohort Calibration Verification (Part 0)
#       - Block 2.1: Stage Hazard Ratios With & Without Full-Cohort LASSO Genes
#   MODULE 3: Unknown-Category Predictive Value (Part 1)
#       - Block 3.1: Paired Bootstrap Delta C-Index & Full-Cohort LRT
#   MODULE 4: Shared-Split Training Comparison (Part 2)
#       - Block 4.1: 100 Shared Stratified Splits on Documented Test Patients
#       - Block 4.2: Nadeau-Bengio Corrected Paired t-Test Summarizer
#   MODULE 5: Execution & Export
#       - Block 5.1: Run Full Pipeline & Save Supplementary CSV Tables (Table S4)
#
# Save as:  12d_inferential_strategy_comparison.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_base.rds
# Output:   TableS_shared_split_A_vs_B.csv, TableS_shared_split_raw.csv
#
# Dependencies: survival, glmnet
#
# REVISION NOTES (Nov 2026):
#   - Refactored into standardized Module -> Block architecture matching Script 12.
#   - Replaced dynamic get()/assign() scoping in Part (1) with explicit data
#     preparation helpers and added a zero-selection guard in Part (0).
#   - Formatted TableS_shared_split_A_vs_B.csv with an explicit 'Model' column
#     (row.names = FALSE) and added execution guard (nsclc.emp12d_execute).
# ==============================================================================


# ==============================================================================
# MODULE 1: Environment & Dependency Management
# ==============================================================================

# --- Block 1.1: Load Required Packages ----------------------------------------
suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
})

# --- Block 1.2: Load Base Cohort & Split Objects ------------------------------
load_base_cohort_objects <- function(base_path = "nsclc_base.rds") {
  if (!file.exists(base_path)) stop("Required input file not found: ", base_path)
  
  base <- readRDS(base_path)
  base$cohort_B_expr$stage_group <- factor(
    base$cohort_B_expr$stage_group,
    levels = c("Early", "Advanced", "Unknown")
  )
  return(base)
}

calc_cindex <- function(time, event, risk) {
  1 - survival::concordance(Surv(time, event) ~ risk)$concordance
}


# ==============================================================================
# MODULE 2: Full-Cohort Calibration Verification (Part 0)
# ==============================================================================

# --- Block 2.1: Stage Hazard Ratios With & Without Full-Cohort LASSO Genes ----
# Note: glmnet standardizes internally; MIR124-1 has zero variance in this cohort,
# so scaling is restricted strictly to non-zero LASSO-selected genes.
run_calibration_check <- function(coh, genes, seed = 123) {
  set.seed(seed)
  x <- as.matrix(coh[, genes])
  y <- Surv(coh$os_time, coh$os_event)
  
  cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, type.measure = "C", nfolds = 10)
  cf     <- as.numeric(coef(cv_fit, s = "lambda.min"))
  sel    <- genes[cf != 0]
  
  cat(
    "\n--- Part (0): Full-Cohort Calibration Model ---\n",
    "Full-cohort LASSO (lambda.min, C-tuned) selected ",
    length(sel), " of ", length(genes), " genes: ",
    paste(sel, collapse = ", "), "\n",
    sep = ""
  )
  
  f_no_genes <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = coh)
  
  if (length(sel) > 0L) {
    z_mat <- as.data.frame(scale(coh[, sel, drop = FALSE]))
    names(z_mat) <- paste0("z_", make.names(sel))
    dd <- cbind(coh, z_mat)
    fmla_with <- as.formula(paste(
      "Surv(os_time, os_event) ~ age + gender + histology + stage_group +",
      paste(names(z_mat), collapse = " + ")
    ))
    f_with_genes <- coxph(fmla_with, data = dd)
  } else {
    f_with_genes <- f_no_genes
  }
  
  extract_stage_hr <- function(fit) {
    round(summary(fit)$conf.int[c("stage_groupAdvanced", "stage_groupUnknown"), c(1, 3, 4)], 2)
  }
  
  cat("Stage HRs without gene adjustment:\n")
  print(extract_stage_hr(f_no_genes))
  cat("Stage HRs with LASSO-selected genes:\n")
  print(extract_stage_hr(f_with_genes))
  
  invisible(list(selected_genes = sel, fit_no_genes = f_no_genes, fit_with_genes = f_with_genes))
}


# ==============================================================================
# MODULE 3: Unknown-Category Predictive Value (Part 1)
# ==============================================================================

# --- Block 3.1: Paired Bootstrap Delta C-Index & Full-Cohort LRT --------------
prepare_binary_stage <- function(df) {
  df$stage_group <- factor(df$stage_group, levels = c("Early", "Advanced", "Unknown"))
  df$stage_bin   <- factor(
    ifelse(df$stage_group == "Advanced", "Advanced", "NotAdvanced"),
    levels = c("NotAdvanced", "Advanced")
  )
  df
}

test_unknown_predictive_value <- function(base, B = 2000, seed = 1) {
  tr  <- prepare_binary_stage(base$train_B_expr)
  te  <- prepare_binary_stage(base$test_B_expr)
  coh <- prepare_binary_stage(base$cohort_B_expr)
  
  m_with_unk <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = tr)
  m_lumped   <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_bin,   data = tr)
  
  r_with <- as.numeric(predict(m_with_unk, newdata = te, type = "lp"))
  r_lump <- as.numeric(predict(m_lumped,   newdata = te, type = "lp"))
  
  c_with <- calc_cindex(te$os_time, te$os_event, r_with)
  c_lump <- calc_cindex(te$os_time, te$os_event, r_lump)
  delta_est <- c_with - c_lump
  
  set.seed(seed)
  n <- nrow(te)
  delta_boot <- replicate(B, {
    idx <- sample.int(n, replace = TRUE)
    calc_cindex(te$os_time[idx], te$os_event[idx], r_with[idx]) -
      calc_cindex(te$os_time[idx], te$os_event[idx], r_lump[idx])
  })
  
  p_boot <- min(1, 2 * min(mean(delta_boot <= 0), mean(delta_boot >= 0)))
  
  cat("\n--- Part (1): Predictive Value of the Unknown Category ---\n")
  cat(
    "Validation C with Unknown category:", round(c_with, 3),
    "| Unknown lumped with Early:", round(c_lump, 3), "\n"
  )
  cat(
    "Delta C =", round(delta_est, 3),
    "95% CI:", paste(round(quantile(delta_boot, c(0.025, 0.975)), 3), collapse = " to "),
    "| p_boot =", signif(p_boot, 3), "\n"
  )
  
  # Full-cohort likelihood-ratio test
  a_with <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = coh)
  a_lump <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_bin,   data = coh)
  lrt <- anova(a_lump, a_with)
  cat("Likelihood-ratio test for Unknown category (full cohort):\n")
  print(lrt)
  
  invisible(list(delta_C = delta_est, delta_boot = delta_boot, lrt = lrt))
}


# ==============================================================================
# MODULE 4: Shared-Split Training Comparison (Part 2)
# ==============================================================================

# --- Block 4.1: 100 Shared Stratified Splits on Documented Test Patients ------
run_shared_splits_comparison <- function(coh, genes, R = 100, base_seed = 1000) {
  out <- matrix(
    NA_real_,
    nrow = R,
    ncol = 4,
    dimnames = list(NULL, c("clin_A", "clin_B", "lasso_A", "lasso_B"))
  )
  
  strata <- interaction(coh$os_event, coh$stage_group, coh$histology, drop = TRUE)
  ntr <- numeric(R)
  nte <- numeric(R)
  
  for (r in seq_len(R)) {
    set.seed(base_seed + r)
    test_idx <- unlist(lapply(split(seq_len(nrow(coh)), strata), function(i) {
      if (length(i) > 1L) sample(i, round(0.3 * length(i))) else integer(0)
    }))
    
    trB <- coh[-test_idx, ]
    trA <- trB[trB$stage_group != "Unknown", ]
    trA$stage_group <- droplevels(trA$stage_group)
    
    doc <- coh[test_idx, ]
    doc <- doc[doc$stage_group != "Unknown", ]
    
    ntr[r] <- nrow(trB)
    nte[r] <- nrow(doc)
    
    for (S in c("A", "B")) {
      trn <- if (S == "A") trA else trB
      
      # Clinical-only Cox model scored on identical documented-stage test patients
      cm <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = trn)
      doc_eval <- doc
      if (S == "A") {
        doc_eval$stage_group <- factor(doc_eval$stage_group, levels = levels(trA$stage_group))
      }
      out[r, paste0("clin_", S)] <- calc_cindex(
        doc$os_time, doc$os_event,
        as.numeric(predict(cm, newdata = doc_eval, type = "lp"))
      )
      
      # Gene-only LASSO Cox model scored on identical documented-stage test patients
      xt  <- as.matrix(trn[, genes])
      cvr <- cv.glmnet(
        x            = xt,
        y            = Surv(trn$os_time, trn$os_event),
        family       = "cox",
        type.measure = "C",
        nfolds       = 10
      )
      pr <- as.numeric(predict(cvr, newx = as.matrix(doc[, genes]), s = "lambda.min"))
      out[r, paste0("lasso_", S)] <- calc_cindex(doc$os_time, doc$os_event, pr)
    }
    
  }
  
  list(out = out, ntr = ntr, nte = nte)
}

# --- Block 4.2: Nadeau-Bengio Corrected Paired t-Test Summarizer --------------
summarize_shared_splits <- function(split_res) {
  out <- split_res$out
  ntr <- split_res$ntr
  nte <- split_res$nte
  R   <- nrow(out)
  rho <- mean(nte / ntr)
  
  nb_paired_test <- function(xB, xA, model_name) {
    d     <- xB - xA
    d_bar <- mean(d)
    s2    <- var(d)
    se_nb <- sqrt((1 / R + rho / (1 - rho)) * s2)
    t_val <- d_bar / se_nb
    p_val <- 2 * pt(-abs(t_val), df = R - 1)
    ci    <- d_bar + c(-1, 1) * qt(0.975, df = R - 1) * se_nb
    
    data.frame(
      Model          = model_name,
      Mean_A         = round(mean(xA), 4),
      Mean_B         = round(mean(xB), 4),
      Mean_Diff_B_A  = round(d_bar, 4),
      NB_SE          = round(se_nb, 4),
      CI_lower       = round(ci[1], 4),
      CI_upper       = round(ci[2], 4),
      t_stat         = round(t_val, 3),
      df             = R - 1,
      p_value        = signif(p_val, 3),
      stringsAsFactors = FALSE
    )
  }
  
  res_tab <- rbind(
    nb_paired_test(out[, "clin_B"],  out[, "clin_A"],  "Clinical Cox"),
    nb_paired_test(out[, "lasso_B"], out[, "lasso_A"], "LASSO Cox")
  )
  
  cat("\n--- Part (2): 100 Shared Splits (Nadeau-Bengio Corrected Paired t-Test) ---\n")
  print(res_tab)
  return(res_tab)
}


# ==============================================================================
# MODULE 5: Execution & Export
# ==============================================================================

# --- Block 5.1: Run Full Pipeline & Save Supplementary CSV Tables (Table S4) --
if (isTRUE(getOption("nsclc.emp12d_execute", TRUE))) {
  cat("========================================================\n")
  cat("Starting Script 12d: Inferential Strategy Comparison\n")
  cat("========================================================\n")
  
  base <- load_base_cohort_objects("nsclc_base.rds")
  
  # Part (0)
  run_calibration_check(base$cohort_B_expr, base$present_genes, seed = 123)
  
  # Part (1)
  test_unknown_predictive_value(base, B = 2000, seed = 1)
  
  # Part (2)
  split_res <- run_shared_splits_comparison(base$cohort_B_expr, base$present_genes, R = 100, base_seed = 1000)
  summary_12d <- summarize_shared_splits(split_res)
  
  write.csv(summary_12d, "TableS_shared_split_A_vs_B.csv", row.names = FALSE)
  write.csv(as.data.frame(split_res$out), "TableS_shared_split_raw.csv", row.names = FALSE)
  
  cat("\nSaved TableS_shared_split_A_vs_B.csv and TableS_shared_split_raw.csv.\n")
  cat("========================================================\n")
  cat("Script 12d Complete!\n")
  cat("========================================================\n")
}
