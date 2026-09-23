source("R/covariance_utils.R")
source("R/simulate.R")
source("R/fit.R")
source("R/harmonize.R")
source("R/standardize.R")
source("R/evaluation_metrics.R")
source("R/ggcombat.R")

set.seed(20260818)

site1 <- simulate_ggcombat_data(
  n = 80,
  block_sizes = c(8, 8, 4),
  latent_cov = matrix(c(1.0, 0.3, 0.3, 1.0), 2, 2),
  batch_means = c(0.4, -0.1, 0.2),
  batch_sds = c(0.15, 0.15, 0.10),
  lambdas = c(8, 8, 8),
  thetas = c(6, 6, 6)
)

site2 <- simulate_ggcombat_data(
  n = 80,
  block_sizes = c(8, 8, 4),
  latent_cov = matrix(c(1.0, -0.2, -0.2, 1.0), 2, 2),
  batch_means = c(-0.3, 0.2, -0.1),
  batch_sds = c(0.15, 0.15, 0.10),
  lambdas = c(10, 7, 9),
  thetas = c(4, 7, 5)
)

Y <- rbind(site1$Z, site2$Z)
batch <- rep(c("site1", "site2"), each = 80)
L <- site1$L

before <- evaluate_site_distances(split(as.data.frame(Y), batch))
before <- evaluate_site_distances(lapply(split(seq_len(nrow(Y)), batch), function(idx) Y[idx, ]))

fit <- ggcombat(Y, batch, L = L, target = "average")

after <- evaluate_site_distances(
  lapply(split(seq_len(nrow(fit$Z_harmonized)), batch), function(idx) fit$Z_harmonized[idx, ])
)

print(before)
print(after)
