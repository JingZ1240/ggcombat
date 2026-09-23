# Revision Robustness Analyses

This folder contains public scripts for reviewer-requested robustness analyses
that can be run without restricted ABCD subject-level data.

Planned analyses:

1. Correlation topology versus covariance magnitude.
2. Operating range under additive and covariance site-effect strength.
3. Block-coherence weakening.
4. Community misspecification.
5. Nonlinear covariate sensitivity.
6. Demographic-site imbalance sensitivity, when expressible using simulated or
   otherwise shareable data.

Implementation rule:

- New scripts should not overwrite manuscript results or private outputs.
- New outputs should be written to clearly named local result folders that are
  excluded from Git tracking.
- Scripts in this public folder should use simulated or otherwise shareable
  inputs.
