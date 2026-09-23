# Nonlinear Covariate Sensitivity

Purpose: address Reviewer 1 Major Comment 2 on nonlinear covariate adjustment
and the relationship between GG-ComBat and DeepComBat.

Planned scripts:

1. `01_run_nonlinear_covariate_sensitivity.R`
   - simulate site-imbalanced age;
   - generate linear, quadratic, and smooth nonlinear biological age effects;
   - compare ComBat, CovBat, and GG-ComBat under linear and spline covariate
     specifications.

2. `02_summarize_nonlinear_covariate_sensitivity.R`
   - summarize site-removal and biological-signal-preservation metrics;
   - create manuscript/supplement candidate figures.

3. `03_deepcombat_feasibility.R`
   - optional feasibility check for DeepComBat installation and a small pilot
     run.

Important scope:

- Keep `L` correctly specified here. Misspecified community structure belongs
  to the separate R2.2 robustness analysis.
- Report biological preservation and site removal separately.
- Do not claim GG-ComBat outperforms DeepComBat unless DeepComBat is actually
  run under a fair, reproducible setup.
