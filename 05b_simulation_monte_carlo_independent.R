# ==============================================================================
# Script 05b: NSCLC Missing Data Simulation Arm - Independent Replicates Extension
#
# Author:   Getachew Derbew Belay
# Date:     November 2026 (revised)
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Extend the Monte Carlo engine (Script 05) with two design changes:
#           (1) Regenerate latent severity AND survival times in every replicate,
#               so replicates are truly independent datasets.
#           (2) Add S0 (full-data reference, no missingness) and S3b
#               (outcome-informed imputation with proper parameter draws).
#           Reuses the missingness generators and evaluate_strategy() from
#           Script 05, and the latent/survival simulators from
#           simulation_helpers.R.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Base Cohort and Shared Code (Script 05 helpers)
#   MODULE 2: S0 Reference and S3b Strategy
#       - Block 2.1: Nelson-Aalen Cumulative Hazard Helper
#       - Block 2.2: S3b Evaluation (proper MI with parameter draws)
#       - Block 2.3: S0 Full-Data Reference
#   MODULE 3: Monte Carlo Engine (Independent Replicates)
#       - Block 3.1: Single-Replicate Execution
#       - Block 3.2: Execution Loop (with Seed)
#   MODULE 4: Summarization & Export
#       - Block 4.1: Calculate Bias, MCSE, RMSE, Coverage
#       - Block 4.2: Save Results
#
# Note: Module 2 of Script 05 (missingness generators gen_mcar/gen_mar/gen_mnar)
#       is reused unchanged; it is not redefined here.
#
# Save as:  05b_simulation_independent_replicates.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_base.rds                (cohort_B_expr)
#           simulation_helpers.R            (assign_latent_stages_param,
#                                            simulate_survival_param)
#           05_simulation_monte_carlo.R     (gen_*, evaluate_strategy;
#                                            sourced with the master gate OFF)
# Output:   nsclc_sim_mc_independent.rds
#           Table4_independent_replicates.csv
#
# Dependencies: survival, dplyr, MASS
#
# REVISION NOTES (Nov 2026):
#   - Replicates are independent: latent severity and survival times are
#     regenerated inside the loop. Script 05 fixed gt and resampled only
#     missingness; this tests the mechanism but underestimates MC variability.
#   - S0 added as an explicit full-data reference fitted to gt before any
#     missingness is introduced. Bias is reported vs. the true value
#     (log 1.86 = 0.620) and vs. S0 (paired within replicate).
#   - S3b added: stochastic imputation using age, histology, event indicator,
#     and Nelson-Aalen H(t), with a parameter draw from the imputation-model
#     posterior (proper-style MI). M = 20 imputations per replicate.
#   - Strategy labels are assigned explicitly inside each evaluator call,
#     not by row position (removes the fragile c("S0","S1",...) binding).
#   - Seeds: 10000+i (latent severity), 20000+i (survival), 30000+i
#     (missingness). Same 30000+i across mechanisms within a replicate, so
#     MCAR/MAR/MNAR differences are attributable to mechanism, not RNG drift.
#   - check: assign_latent_stages_param(unk_hr = 3.22, ...) assumes the helper
#     takes an HR and applies log() internally. If the helper expects a log
#     hazard, change to unk_hr = log(3.22). Verify in simulation_helpers.R
#     before the production run; Sections 2.5 and 2.6 of the manuscript
#     should use the same convention.
#
# USAGE:
#   Rscript 05b_simulation_independent_replicates.R
#   Smoke test first: run_simulation_pipeline(n_iter = 5).
#   Outputs are written to the current working directory (same as Script 05).
#   Run from the NSCLC_Missing_Data folder so outputs land beside the script.
# ==============================================================================

# ##############################################################################
# MODULE 1: Environment & Dependency Management
# ##############################################################################

# --- Block 1.1: Load Required Packages ---
init_environment_05b <- function() {
  cat("\n>>> MODULE 1, BLOCK 1.1: Loading Packages\n")
  
  packages <- c("survival", "dplyr", "MASS")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# --- Block 1.2: Load Base Cohort and Shared Code (Script 05 helpers) ---
load_base_cohort <- function() {
  cat(">>> MODULE 1, BLOCK 1.2: Loading Base Cohort and Shared Code\n")
  
  if (!file.exists("nsclc_base.rds")) {
    stop("Run the NSCLC Missing Data preparation script first to generate nsclc_base.rds.")
  }
  
  base   <- readRDS("nsclc_base.rds")
  cohort <- base$cohort_B_expr
  
  # Source shared helpers and Script 05 (generators + evaluate_strategy).
  # 05's master gate is temporarily switched OFF so that sourcing 05 does not
  # re-run its full Monte Carlo. The caller's original option is restored.
  old_opt <- getOption("nsclc.mc_execute", default = TRUE)
  
  if (!file.exists("simulation_helpers.R")) {
    stop("simulation_helpers.R not found in working directory.")
  }
  source("simulation_helpers.R")
  
  if (!file.exists("05_simulation_monte_carlo.R")) {
    stop("05_simulation_monte_carlo.R not found in working directory.")
  }
  options(nsclc.mc_execute = FALSE)
  source("05_simulation_monte_carlo.R")
  options(nsclc.mc_execute = old_opt)
  
  true_beta_adv <- log(1.86)
  
  return(list(cohort        = cohort,
              true_beta_adv = true_beta_adv))
}


# ##############################################################################
# MODULE 2: S0 Reference and S3b Strategy
# ##############################################################################

# --- Block 2.1: Nelson-Aalen Cumulative Hazard Helper ---
# H(t) is computed on the full simulated data and used as a proxy for
# unobserved severity in the S3b imputation model. method = "constant" and
# rule = 2 give a right-continuous step function that extrapolates flat beyond
# the last observed event time.
nelsonaalen_H <- function(time, event) {
  fit <- survfit(Surv(time, event) ~ 1)
  approx(fit$time, fit$cumhaz,
         xout = time, method = "constant", f = 0, rule = 2)$y
}

# --- Block 2.2: S3b Evaluation (Outcome-Informed Proper MI) ---
# Imputation model: logistic regression of Advanced vs Early on age, histology,
# event indicator, and Nelson-Aalen H(t).
#
# For each of M imputations:
#   (1) draw a parameter vector from N(hat_beta, V_hat_beta)
#   (2) draw imputed stage labels from the resulting predicted probabilities
#   (3) fit Cox on the completed data
#
# Pooled estimate uses Rubin's rules: mean estimate; SE combines within- and
# between-imputation variance with the (1 + 1/M) factor.
#
# Including the outcome in the imputation model is not leakage: the imputation
# model is fit on the observed-stage subset only, and the event indicator and
# H(t) are observed for all subjects at analysis time. Omitting the outcome
# attenuates Cox hazard ratios toward the null (White & Royston 2009,
# Stat Med 28:1982).
evaluate_S3b <- function(data, M = 20) {
  tryCatch({
    required <- c("sim_time", "sim_event", "age", "histology", "observed_stage")
    missing_cols <- setdiff(required, names(data))
    if (length(missing_cols) > 0) {
      stop("evaluate_S3b: missing columns: ", paste(missing_cols, collapse = ", "))
    }
    
    d <- data[, required]
    d$H <- nelsonaalen_H(d$sim_time, d$sim_event)
    
    obs <- d[d$observed_stage != "Unknown", ]
    mis <- d[d$observed_stage == "Unknown", ]
    
    if (nrow(mis) == 0) {
      return(data.frame(beta_adv = NA, se_adv = NA, beta_unknown = NA))
    }
    
    obs$adv <- as.integer(obs$observed_stage == "Advanced")
    fm <- glm(adv ~ age + histology + sim_event + H,
              data = obs, family = binomial)
    
    beta_adv_imps <- numeric(M)
    se_adv_imps   <- numeric(M)
    
    for (m in seq_len(M)) {
      th <- MASS::mvrnorm(1, coef(fm), vcov(fm))
      X  <- model.matrix(~ age + histology + sim_event + H, data = mis)
      p  <- plogis(drop(X %*% th))
      
      comp <- rbind(obs[, names(mis)], mis)
      comp$observed_stage <- as.character(comp$observed_stage)
      comp$observed_stage[comp$observed_stage == "Unknown"] <-
        ifelse(runif(nrow(mis)) < p, "Advanced", "Early")
      comp$observed_stage <- factor(comp$observed_stage,
                                    levels = c("Early", "Advanced"))
      
      fit <- coxph(Surv(sim_time, sim_event) ~ observed_stage, data = comp)
      beta_adv_imps[m] <- coef(fit)["observed_stageAdvanced"]
      se_adv_imps[m]   <- sqrt(vcov(fit)["observed_stageAdvanced", "observed_stageAdvanced"])
    }
    
    # Rubin's rules
    beta_adv <- mean(beta_adv_imps)
    se_adv   <- sqrt(mean(se_adv_imps^2) + (1 + 1/M) * var(beta_adv_imps))
    
    return(data.frame(
      beta_adv     = as.numeric(beta_adv),
      se_adv       = as.numeric(se_adv),
      beta_unknown = NA_real_
    ))
  }, error = function(e) {
    return(data.frame(beta_adv = NA_real_, se_adv = NA_real_, beta_unknown = NA_real_))
  })
}

# --- Block 2.3: S0 Full-Data Reference ---
# Fitted to gt before any missingness is introduced. The estimand is the
# Advanced-vs-Early contrast among documented patients; the 123 original
# Unknown patients form a third latent-risk level ("Ultra-High").
evaluate_S0 <- function(gt) {
  tryCatch({
    fit <- coxph(Surv(sim_time, sim_event) ~ stage_group, data = gt)
    return(data.frame(
      beta_adv     = as.numeric(coef(fit)["stage_groupAdvanced"]),
      se_adv       = as.numeric(sqrt(vcov(fit)["stage_groupAdvanced", "stage_groupAdvanced"])),
      beta_unknown = NA_real_
    ))
  }, error = function(e) {
    return(data.frame(beta_adv = NA_real_, se_adv = NA_real_, beta_unknown = NA_real_))
  })
}


# ##############################################################################
# MODULE 3: Monte Carlo Engine (Independent Replicates)
# ##############################################################################

# --- Block 3.1: Single-Replicate Execution ---
run_one_replicate <- function(cohort, mech, i) {
  # (a) regenerate latent severity and survival times
  gt <- cohort |>
    assign_latent_stages_param(unk_hr = 3.22, seed = 10000 + i) |>
    simulate_survival_param(adv_hr = 1.86, early_median = 1984, seed = 20000 + i)
  
  # (b) inject missingness under the current mechanism
  set.seed(30000 + i)
  sim_data <- switch(mech,
                     MCAR = gen_mcar(gt),
                     MAR  = gen_mar(gt),
                     MNAR = gen_mnar(gt))
  
  # (c) apply all strategies to the same dataset
  rows <- list()
  
  res <- evaluate_S0(gt)
  rows[[length(rows) + 1]] <- data.frame(
    Mechanism = mech, Iteration = i, Strategy = "S0",
    Beta_Adv = res$beta_adv, SE_Adv = res$se_adv, Beta_Unknown = res$beta_unknown
  )
  
  for (strat in c("S1", "S2", "S3", "S4")) {
    res <- evaluate_strategy(sim_data, strat)
    rows[[length(rows) + 1]] <- data.frame(
      Mechanism = mech, Iteration = i, Strategy = strat,
      Beta_Adv = res$beta_adv, SE_Adv = res$se_adv, Beta_Unknown = res$beta_unknown
    )
  }
  
  res <- evaluate_S3b(sim_data)
  rows[[length(rows) + 1]] <- data.frame(
    Mechanism = mech, Iteration = i, Strategy = "S3b",
    Beta_Adv = res$beta_adv, SE_Adv = res$se_adv, Beta_Unknown = res$beta_unknown
  )
  
  do.call(rbind, rows)
}

# --- Block 3.2: Execution Loop ---
run_monte_carlo_independent <- function(cohort, n_iter = 100) {
  cat("\n>>> MODULE 3, BLOCK 3.2: Running Independent-Replicate Monte Carlo\n")
  cat(sprintf("Iterations per mechanism: %d\n", n_iter))
  
  results <- list()
  
  for (mech in c("MCAR", "MAR", "MNAR")) {
    cat("\n  - Simulating mechanism:", mech, "\n")
    
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    
    for (i in 1:n_iter) {
      setTxtProgressBar(pb, i)
      results[[length(results) + 1]] <- run_one_replicate(cohort, mech, i)
    }
    close(pb)
  }
  
  results_df <- do.call(rbind, results)
  return(results_df)
}


# ##############################################################################
# MODULE 4: Summarization & Export
# ##############################################################################

# --- Block 4.1: Calculate Bias, MCSE, RMSE, Coverage ---
summarize_mc_results_independent <- function(results_df, true_beta_adv) {
  cat("\n>>> MODULE 4, BLOCK 4.1: Summarizing Results\n")
  
  # Paired S0 reference per (Mechanism, Iteration)
  s0 <- results_df[results_df$Strategy == "S0",
                   c("Mechanism", "Iteration", "Beta_Adv")]
  names(s0)[3] <- "Beta_S0"
  
  merged <- merge(results_df, s0, by = c("Mechanism", "Iteration"))
  
  summary_df <- merged %>%
    group_by(Mechanism, Strategy) %>%
    summarize(
      Bias_Adv        = mean(Beta_Adv - true_beta_adv, na.rm = TRUE),
      MCSE_bias       = sd(Beta_Adv, na.rm = TRUE) / sqrt(sum(!is.na(Beta_Adv))),
      Bias_vs_S0      = mean(Beta_Adv - Beta_S0, na.rm = TRUE),
      Empirical_SD    = sd(Beta_Adv, na.rm = TRUE),
      Mean_model_SE   = mean(SE_Adv, na.rm = TRUE),
      RMSE            = sqrt(mean((Beta_Adv - true_beta_adv)^2, na.rm = TRUE)),
      Coverage        = mean((Beta_Adv - 1.96 * SE_Adv <= true_beta_adv) &
                               (Beta_Adv + 1.96 * SE_Adv >= true_beta_adv), na.rm = TRUE),
      .groups = "drop"
    )
  
  print(as.data.frame(summary_df))
  return(summary_df)
}

# --- Block 4.2: Save Results ---
# Writes both the RDS and the CSV table to the current working directory,
# matching Script 05's convention. No figure is produced by this script;
# the bias/coverage figure is generated by the plotting script.
save_results_independent <- function(raw_results, summary_df) {
  cat("\n>>> MODULE 4, BLOCK 4.2: Saving Results\n")
  
  saveRDS(list(
    raw_results = raw_results,
    summary     = summary_df
  ), "nsclc_sim_mc_independent.rds")
  
  write.csv(summary_df, "Table4_independent_replicates.csv", row.names = FALSE)
  
  cat("Saved nsclc_sim_mc_independent.rds successfully.\n")
  cat("Saved Table4_independent_replicates.csv successfully.\n")
}


# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_simulation_pipeline <- function(n_iter = 100) {
  cat("========================================================\n")
  cat("Starting NSCLC Missing Data Independent-Replicate Monte Carlo\n")
  cat("========================================================\n")
  
  init_environment_05b()
  setup <- load_base_cohort()
  
  raw_results <- run_monte_carlo_independent(setup$cohort, n_iter = n_iter)
  
  summary_df <- summarize_mc_results_independent(raw_results, setup$true_beta_adv)
  
  save_results_independent(raw_results, summary_df)
  
  cat("\n========================================================\n")
  cat("Simulation Pipeline Complete!\n")
  cat("========================================================\n")
}

# Execute unless a sourcing script has opted out. Scripts that source 05b
# only for its helpers should set options(nsclc.mc_execute = FALSE).
if (isTRUE(getOption("nsclc.mc_execute", default = TRUE))) {
  run_simulation_pipeline()
}
