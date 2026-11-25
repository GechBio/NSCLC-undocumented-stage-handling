# Handling undocumented pathologic stage in NSCLC prognostic modelling: analysis code

This repository contains the R pipeline that reproduces the results in the manuscript
*"Handling Undocumented Pathologic Stage in NSCLC Prognostic Models Built from Public Data:
Complete-Case Exclusion versus Retention in a TCGA Analysis and Simulation Study"*

## Contents

| File | Purpose |
|---|---|
| `00_base_setup.R` | Cohort assembly from TCGA-LUAD/LUSC (N = 867 expression, N = 701 methylation), train/validation splits |
| `01_empirical_expression.R`, `02_empirical_methylation.R` | Empirical models, Strategy A (complete-case) vs B (Unknown retained) |
| `03_generate_empirical_figures.R`, `03b_feature_selection_upset.R` | Figures 2, 3, and 4 |
| `04_simulation_ground_truth.R`, `04b_generate_simulation_ground_truth_figure.R` | Calibrated simulation ground truth (Table 3, Figure 5) |
| `05_simulation_monte_carlo.R` | Monte Carlo engine (functions sourced by 05b and 05c; not a pipeline step) |
| `05b_simulation_monte_carlo_independent.R` | 100 independent replicates under MCAR, MAR, MNAR (Table 4) |
| `06b_script_additional.R` | Figure 6 |
| `05c_simulation_stage_dependent_mnar.R` | Four stage-dependent MNAR scenarios (Table S5) |
| `07_reporting_checklist.R` | Table 5 |
| `09_clinical_cohort_disambiguation.R` | Tables S1 and S2 (Kaplan-Meier median vs observation time and consistency checks) |
| `script12_additional.R`, `12b_script_additional.R`, `script_12c.R`, `Script_12d_Inferential_.R`, `12e_rsf_pvalues.R` | Bootstrap C-index CIs, 100 repeated splits and selection frequency, paired and shared-split tests, RSF p-values (Tables S3–S4) |
| `simulation_helpers.R` | Helper functions for the simulation |
| `08_master_pipeline.R` | Runs everything in the correct order |
| `renv.lock` | Frozen R package environment |

## Folder layout (important)

Script `00_base_setup.R` expects this layout and writes outputs inside `NSCLC_Missing_Data/`:

```text
<project root>/
├── NSCLC_GDC_Sample_Sheet.tsv      (GDC sample sheet)
├── master_nsclc_matrix.rds         (expression matrix)
├── clinical_final.rds              (survival endpoints)
├── clinical_covariates.rds         (age, gender, stage, project)
├── nsclc_raw_data/                 (*.methylation_array.sesame.level3betas.txt)
└── NSCLC_Missing_Data/                       <- this repository (via git clone or ZIP download)
```

Methylation beta files (`*.methylation_array.sesame.level3betas.txt`) go in `<project root>/nsclc_raw_data`
or in a folder given by the environment variable `METH_DIR`. `PROJECT_ROOT` can be set to change the root.
To skip the raw-data steps, place the provided `nsclc_base.rds` in `NSCLC_Missing_Data/` and comment out
`"00_base_setup.R"` in `08_master_pipeline.R`.

## Environment

1. Install R 4.5.2.
2. In R, from the repository folder:
   ```r
   install.packages("renv")
   renv::restore()
   ```

## Data

Raw data are not included. Download RNA-seq and HumanMethylation450K data for TCGA-LUAD and TCGA-LUSC
from the NCI Genomic Data Commons (https://cancer.gov) and place them as described in
the header of `00_base_setup.R`. Processed objects (`nsclc_*.rds`) are provided in the
Zenodo archive so that downstream scripts can be run without repeating the download.

## Running

In R, set the working directory to `NSCLC_Missing_Data/` and run:

```r
source("08_master_pipeline.R")
```

The simulation scripts (`05b`, `05c`) take several minutes to hours depending on hardware. Seeds are set
inside each script, so results are reproducible.
