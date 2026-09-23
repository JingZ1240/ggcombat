# Common Study Helpers

This folder contains cleaned helper scripts extracted from the submitted notebooks.

These helpers are intended to preserve the submitted analyses while making the workflow easier to rerun. They should not silently change estimands, default parameters, or output definitions.

Important legacy detail:

- The submitted main simulation uses `G_vec = c(50, 50, 50, 50)` but some GG-ComBat code represents the first three communities explicitly and treats the final 50 features as singleton/weakly connected by position.
- Do not change this representation in reproduction scripts unless the regenerated results are intentionally marked as revised.
