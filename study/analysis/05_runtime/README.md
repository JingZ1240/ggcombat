# Runtime Benchmark

Role in revision:

- Respond to the reviewer request for computational complexity and practical runtime.

Timing should report:

1. Community/graph construction for `L`, separately.
2. Standardization.
3. GG-ComBat empirical-Bayes fitting conditional on fixed `L`.
4. GG-ComBat harmonization conditional on fixed `L`.
5. ComBat total time.
6. CovBat total time.
7. Block-ComBat total time.

Reason for separate `L` timing:

- `L` can be estimated once and reused, so graph construction should not be collapsed into every harmonization run.

Target dimensions:

- Simulation scale: `N = 900`, `G = 200`, `I = 3`.
- ABCD scale: `N = 6467`, `G = 426`, `I = 12`.
