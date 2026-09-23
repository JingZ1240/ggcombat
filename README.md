# ggcombat

`ggcombat` implements Graph-Guided ComBat (GG-ComBat), a statistical
harmonization method for high-dimensional multi-site data with structured
cross-feature dependence. The method was developed for regional neuroimaging
features, but the model is applicable to other biomedical data types in which
features form correlated modules, pathways, subnetworks, or other interpretable
groups.

GG-ComBat extends the empirical-Bayes formulation of ComBat by allowing site
effects to borrow strength within graph-derived feature communities and by
aligning structured residual covariance across sites. In contrast to purely
featurewise harmonization, the method targets both marginal site effects and
second-order site effects when the covariance structure is scientifically
meaningful.

## Statistical Setting

Let `Y` be an `n x p` feature matrix, `batch` a site/scanner indicator, and
`X` optional biological or demographic covariates to preserve during
harmonization. GG-ComBat uses a supplied feature-membership matrix `L` to define
groups of correlated features. Conditional on `L`, the method estimates:

- site-specific additive effects;
- site-specific marginal scale effects;
- block-guided empirical-Bayes hyperparameters;
- site-specific structured covariance components.

The working assumption is not that all sites have identical covariance before
harmonization. Rather, GG-ComBat assumes that a common feature grouping is a
useful scaffold for borrowing information, while allowing site-specific means,
marginal variances, and within-block covariance parameters.

Singleton or weakly connected features can be handled by featurewise
ComBat-style updates, so the method does not force all variables into
communities.

## Installation

```r
install.packages("devtools")
devtools::install_local("ggcombat")
```

During development, from the package root:

```r
devtools::load_all(".")
devtools::test()
```

## Minimal Example

```r
source("examples/quickstart_simulation.R")
```

The example simulates two sites with block-structured covariance, fits
GG-ComBat, and compares site-pair distances before and after harmonization.

## Public API

```r
fit <- ggcombat(
  Y,
  batch,
  covariates = NULL,
  L,
  target = "average"
)
```

Main arguments:

- `Y`: numeric subject-by-feature matrix.
- `batch`: site, scanner, study, or other batch indicator.
- `covariates`: optional data frame of biological covariates to preserve.
- `L`: binary or working membership matrix defining graph-guided feature groups.
- `target`: covariance target, currently including an average target.

The fitted object contains the harmonized matrix and intermediate quantities
used in standardization, empirical-Bayes estimation, and covariance alignment.

## Harmonization Map

After covariate standardization and additive site-effect removal, the
correlated-feature transformation is based on the decomposition

```text
Sigma_i = D_i R_i D_i
```

where `D_i` is diagonal and `R_i` is the site-specific correlation matrix within
the modeled feature block. The row-vector transformation is:

```text
(Z_i - gamma_i) D_i^{-1} R_i^{-1/2} R_target^{1/2}
```

The submitted implementation then applies empirical marginal calibration to keep
the realized post-alignment marginal variances close to the target in finite
samples. For sensitivity analyses, `variance_calibration = "model_target"` uses
the direct model-based target scale. Singleton features use univariate
ComBat-style standardization by default.

## Choosing the Graph Scaffold

GG-ComBat conditions on the supplied matrix `L`. In applications, `L` may come
from anatomical groupings, known networks, pathway annotations, or data-adaptive
community detection. The manuscript analyses use dense subgraph discovery on a
weighted correlation matrix computed from standardized residual features. This
choice is intended to estimate correlation topology after removing covariate and
site-location components, rather than using raw feature covariance magnitudes.

Because uncertainty in `L` is not propagated by the current estimator, applied
analyses should report sensitivity to reasonable graph-construction choices.

## Repository Map

- `R/`: reusable GG-ComBat functions.
- `man/`: help files for the public API.
- `examples/`: small synthetic examples.
- `tests/testthat/`: package-level regression tests.
- `docs/method_notes.md`: implementation notes linking code to the method.
- `study/analysis/01_main_simulation/`: simulation code for the main paper setting.
- `study/analysis/02_np_parameter_recovery/`: simulation code for parameter-recovery settings.
- `study/analysis/04_revision_robustness/`: reviewer-requested robustness simulation scripts.
- `paper/abcd/README.md`: restricted ABCD input schema and data-use notes.

The public repository is intentionally package-centered. It includes the method
implementation, tests, toy examples, and simulation scripts needed to understand
and reproduce the numerical settings described in the manuscript. It does not
include raw ABCD data, subject-level derived ABCD files, private review
materials, manuscript source files, or exploratory notebooks.

## Reproducibility And Data Access

The repository includes reusable method code, simulation scripts, package tests,
and synthetic example data. Raw ABCD data and executable ABCD application
workflows are not included because the ABCD data are restricted by data-use
agreements and the subject-level preprocessing environment is not public.
`paper/abcd/README.md` documents the data-access boundary and the expected
private input types for authorized users.

Do not commit raw ABCD files, subject-level derived files, large `.RData`
workspaces, local HTML renders, private review materials, manuscript drafts, or
credentials.

## Current Status

The package build and tests are checked locally with:

```r
devtools::test()
```

and:

```bash
R CMD build .
R CMD check --no-manual --no-vignettes ggcombat_0.0.1.tar.gz
```

The current public code is intended as a transparent research implementation for
the GG-ComBat manuscript and its revision analyses.
