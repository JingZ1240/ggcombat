# Analysis Scripts

This folder contains public analysis scripts that do not require restricted
subject-level data.

## Public Components

1. `01_main_simulation/`
   - simulation code for the main paper setting;
   - formal entry points are `.R` scripts, not archived notebooks.

2. `02_np_parameter_recovery/`
   - parameter-recovery and varying-dimension simulation code.

3. `04_revision_robustness/`
   - reviewer-requested robustness simulations, including operating-range and
     nonlinear-covariate sensitivity analyses.

4. `05_runtime/`
   - runtime-analysis notes.

## Restricted Components

The ABCD application uses restricted subject-level data and private derived
files. Those workflows are not distributed in this public repository. See
`paper/abcd/README.md` for the data-access boundary.

Generated results, figures, private notebooks, and old exploratory scripts are
excluded from Git tracking.
