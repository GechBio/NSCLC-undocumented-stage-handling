# ==============================================================================
# Script 05c: NSCLC Missing Data Simulation Arm - Stage-Dependent MNAR Boundary Analysis
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Benchmark missing-data handling strategies (S1-S4, S3b) against the
#           full-data reference (S0) under stage-dependent MNAR mechanisms,
#           testing the boundary where complete-case analysis (S1) breaks down
#           when missingness depends jointly on true stage and outcome.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages & Helper Scripts
#       - Block 1.2: Load Base Cohort Data
#   MODULE 2: Stage-Dependent MNAR Generators (Two-Component Design)
#       - Block 2.1: Scenario-Specific Linear Predictor & Calibration
#   MODULE 3: Handling Strategies (S0, S1, S2, S3, S3b, S4)
#       - Block 3.1: Nelson-Aalen Cumulative Hazard Helper
#       - Block 3.2: Outcome-Informed Multiple Imputation (S3b)
#   MODULE 4: Monte Carlo Engine
#       - Block 4.1: Independent Replicate Execution Loop (n_iter = 100)
#   MODULE 5: Summarization & Export
#       - Block 5.1: Calculate Bias, MCSE, SD, RMSE, and 95% CI Coverage
#       - Block 5.2: Save RDS Object & Supplementary Table S5 CSV
#
# Save as:  05c_simulation_stage_mnar.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/NSCLC_Missing_Data/
#
# Input:    nsclc_base.rds, simulation_helpers.R, 05_simulation_monte_carlo.R
# Output:   nsclc_sim_stage_mnar.rds, Table_S5_stage_dependent_MNAR.csv
#
# Dependencies: survival, dplyr, MASS
#
# REVISION NOTES (Nov 2026):
#   - Refactored into standardized Module -> Block architecture matching Script 05.
#   - Evaluates four stage-dependent MNAR boundary scenarios:
#       (A) A_stage_only: missingness depends on Advanced stage (+1.0 * adv)
#       (B) B_advanced_x_death: moderate stage x event interaction (+1.5 * adv * ev)
#       (C) C_early_x_death: early-stage x event interaction (+1.5 * ear * ev)
#       (D) D_advanced_x_death_strong: strong stage x event interaction (+3.0 * adv * ev)
#   - Uses independent replicate survival regeneration (100 replicates/scenario)
#     and aligns strategy output order (S0, S1, S2, S3, S3b, S4) with Table 4.
# ==============================================================================


# ==============================================================================
# MODULE 1: Environment & Dependency Management
# ==============================================================================

# --- Block 1.1: Load Required Packages & Helper Scripts -----------------------
suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(MASS)
})

source("simulation_helpers.R")

# Source Script 05 without triggering its execution loop to reuse evaluate_strategy()
options(nsclc.mc_execute = FALSE)
source("05_simulation_monte_carlo.R")
options(nsclc.mc_execute = TRUE)

# --- Block 1.2: Load Base Cohort Data -----------------------------------------
load_base_cohort <- function(rds_path = "nsclc_base.rds") {
  if (!file.exists(rds_path)) {
    stop("Required input file not found: ", rds_path)
  }
  base <- readRDS(rds_path)
  return(base$cohort_B_expr)
}


# ==============================================================================
# MODULE 2: Stage-Dependent MNAR Generators (Two-Component Design)
# ==============================================================================

# --- Block 2.1: Scenario-Specific Linear Predictor & Calibration --------------
# Model: logit P(missing) = g0 - 0.5 * log(time + 1) + latent_severity + extra
# Intercept g0 is solved via uniroot() so the injected marginal missingness rate
# equals target_rate (14.2%, matching the empirical TCGA cohort and Script 05b).
gen_mnar_stage <- function(data, scenario, target_rate = 0.142) {
  adv <- as.numeric(data$stage_group == "Advanced")
  ear <- as.numeric(data$stage_group == "Early")
  ev  <- as.numeric(data$sim_event)
  
  extra <- switch(
    scenario,
    "A_stage_only"              = 1.0 * adv,       # Advanced stage more often undocumented
    "B_advanced_x_death"        = 1.5 * adv * ev,  # Moderate Advanced x death interaction
    "C_early_x_death"           = 1.5 * ear * ev,  # Moderate Early x death interaction
    "D_advanced_x_death_strong" = 3.0 * adv * ev,  # Strong Advanced x death interaction
    stop("Unknown stage-dependent MNAR scenario: ", scenario)
  )
  
  base_lp <- -0.5 * log(data$sim_time + 1) + data$latent_severity + extra
  g0 <- uniroot(
    function(g) mean(plogis(g + base_lp)) - target_rate,
    interval = c(-60, 20)
  )$root
  
  prob_miss <- plogis(g0 + base_lp)
  data$observed_stage <- as.character(data$stage_group)
  data$observed_stage[runif(nrow(data)) < prob_miss] <- "Unknown"
  data$observed_stage <- factor(
    data$observed_stage,
    levels = c("Early", "Advanced", "Unknown")
  )
  
  return(data)
}


# ==============================================================================
# MODULE 3: Handling Strategies (S0, S1, S2, S3, S3b, S4)
# ==============================================================================

# --- Block 3.1: Nelson-Aalen Cumulative Hazard Helper -------------------------
nelsonaalen_H <- function(time, event) {
  fit <- survfit(Surv(time, event) ~ 1)
  approx(
    x      = fit$time,
    y      = fit$cumhaz,
    xout   = time,
    method = "constant",
    f      = 0,
    rule   = 2
  )$y
}

# --- Block 3.2: Outcome-Informed Multiple Imputation (S3b) --------------------
# Imputes missing stage using age, histology, event indicator, and Nelson-Aalen
# cumulative hazard H(t), drawing parameters from their approximate posterior
# (proper-style MI) and pooling across M = 20 imputations via Rubin's rules.
evaluate_S3b <- function(data, M = 20) {
  tryCatch({
    d <- data[, c("sim_time", "sim_event", "age", "histology", "observed_stage")]
    d$H <- nelsonaalen_H(d$sim_time, d$sim_event)
    
    obs <- d[d$observed_stage != "Unknown", ]
    mis <- d[d$observed_stage == "Unknown", ]
    obs$adv <- as.integer(obs$observed_stage == "Advanced")
    
    fm <- glm(adv ~ age + histology + sim_event + H, data = obs, family = binomial)
    b  <- numeric(M)
    v  <- numeric(M)
    
    for (m in seq_len(M)) {
      th <- MASS::mvrnorm(1, coef(fm), vcov(fm))
      X  <- model.matrix(~ age + histology + sim_event + H, data = mis)
      p  <- plogis(X %*% th)
      
      comp <- rbind(obs[, names(mis)], mis)
      comp$observed_stage <- as.character(comp$observed_stage)
      comp$observed_stage[comp$observed_stage == "Unknown"] <-
        ifelse(runif(nrow(mis)) < p, "Advanced", "Early")
      comp$observed_stage <- factor(comp$observed_stage, levels = c("Early", "Advanced"))
      
      f    <- coxph(Surv(sim_time, sim_event) ~ observed_stage, data = comp)
      b[m] <- coef(f)["observed_stageAdvanced"]
      v[m] <- vcov(f)["observed_stageAdvanced", "observed_stageAdvanced"]
    }
    
    pooled_beta <- mean(b)
    pooled_se   <- sqrt(mean(v) + (1 + 1 / M) * var(b))
    
    data.frame(
      beta_adv     = pooled_beta,
      se_adv       = pooled_se,
      beta_unknown = NA_real_
    )
  }, error = function(e) {
    data.frame(beta_adv = NA_real_, se_adv = NA_real_, beta_unknown = NA_real_)
  })
}


# ==============================================================================
# MODULE 4: Monte Carlo Engine
# ==============================================================================

# --- Block 4.1: Independent Replicate Execution Loop --------------------------
run_stage_mnar <- function(cohort, n_iter = 100) {
  scenarios <- c(
    "A_stage_only",
    "B_advanced_x_death",
    "C_early_x_death",
    "D_advanced_x_death_strong"
  )
  strategy_labels <- c(
    "S0: Full data (reference)",
    "S1: Complete-case",
    "S2: Unknown retention",
    "S3: Imputation, no outcome",
    "S3b: Imputation, + outcome",
    "S4: IPW"
  )
  
  out <- vector("list", length(scenarios) * n_iter)
  idx <- 1L
  
  for (scen in scenarios) {
    message(sprintf("[Script 05c] Running scenario: %s (%d replicates)...", scen, n_iter))
    
    for (i in seq_len(n_iter)) {
      gt <- cohort |>
        assign_latent_stages_param(unk_hr = 3.22, seed = 10000 + i) |>
        simulate_survival_param(adv_hr = 1.86, early_median = 1984, seed = 20000 + i)
      
      set.seed(30000 + i)
      d <- gen_mnar_stage(gt, scen)
      
      # S0: Full-data reference
      f0 <- coxph(Surv(sim_time, sim_event) ~ stage_group, data = gt)
      res_S0 <- data.frame(
        beta_adv     = coef(f0)["stage_groupAdvanced"],
        se_adv       = sqrt(vcov(f0)["stage_groupAdvanced", "stage_groupAdvanced"]),
        beta_unknown = NA_real_
      )
      
      # S1-S4 and S3b
      res_S1  <- evaluate_strategy(d, "S1")
      res_S2  <- evaluate_strategy(d, "S2")
      res_S3  <- evaluate_strategy(d, "S3")
      res_S3b <- evaluate_S3b(d)
      res_S4  <- evaluate_strategy(d, "S4")
      
      r <- rbind(res_S0, res_S1, res_S2, res_S3, res_S3b, res_S4)
      r$Strategy  <- strategy_labels
      r$Mechanism <- scen
      r$Iteration <- i
      
      out[[idx]] <- r
      idx <- idx + 1L
    }
  }
  
  do.call(rbind, out)
}


# ==============================================================================
# MODULE 5: Summarization & Export
# ==============================================================================

# --- Block 5.1: Calculate Bias, MCSE, SD, RMSE, and 95% CI Coverage -----------
summarize_indep <- function(raw) {
  true_beta <- log(1.86)
  
  s0 <- raw[raw$Strategy == "S0: Full data (reference)", c("Mechanism", "Iteration", "beta_adv")]
  names(s0)[3] <- "beta_S0"
  raw <- merge(raw, s0, by = c("Mechanism", "Iteration"))
  
  summary_df <- raw |>
    group_by(Mechanism, Strategy) |>
    summarize(
      Bias_vs_true  = mean(beta_adv - true_beta, na.rm = TRUE),
      MCSE_bias     = sd(beta_adv, na.rm = TRUE) / sqrt(sum(!is.na(beta_adv))),
      Bias_vs_S0    = mean(beta_adv - beta_S0, na.rm = TRUE),
      Empirical_SD  = sd(beta_adv, na.rm = TRUE),
      Mean_model_SE = mean(se_adv, na.rm = TRUE),
      RMSE          = sqrt(mean((beta_adv - true_beta)^2, na.rm = TRUE)),
      Coverage      = mean((beta_adv - 1.96 * se_adv <= true_beta) &
                             (beta_adv + 1.96 * se_adv >= true_beta), na.rm = TRUE),
      .groups       = "drop"
    )
  
  print(as.data.frame(summary_df))
  return(summary_df)
}

# --- Block 5.2: Save RDS Object & Supplementary Table S5 CSV --------------------
if (isTRUE(getOption("nsclc.mnar05c_execute", TRUE))) {
  cat("========================================================\n")
  cat("Starting Script 05c: Stage-Dependent MNAR Analysis\n")
  cat("========================================================\n")
  
  if (!exists("cohort")) {
    base_obj <- readRDS("nsclc_base.rds")
    cohort   <- base_obj$cohort_B_expr
  }
  
  raw_mnar     <- run_stage_mnar(cohort, n_iter = 100)
  summary_mnar <- summarize_indep(raw_mnar)
  
  saveRDS(list(
    raw_results = raw_mnar,
    summary     = summary_mnar
  ), "nsclc_sim_stage_mnar.rds")
  
  write.csv(summary_mnar, "Table_S5_stage_dependent_MNAR.csv", row.names = FALSE)
  
  cat("\nSaved nsclc_sim_stage_mnar.rds and Table_S5_stage_dependent_MNAR.csv.\n")
  cat("========================================================\n")
  cat("Script 05c Complete!\n")
  cat("========================================================\n")
}
