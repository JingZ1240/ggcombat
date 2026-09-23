test_that("ggcombat returns matrices with original dimensions", {
  set.seed(1)
  site1 <- simulate_ggcombat_data(
    n = 20,
    block_sizes = c(4, 4, 2),
    latent_cov = diag(2),
    batch_means = c(0.2, 0.1, 0),
    batch_sds = rep(0.1, 3),
    lambdas = rep(8, 3),
    thetas = rep(6, 3)
  )
  site2 <- simulate_ggcombat_data(
    n = 20,
    block_sizes = c(4, 4, 2),
    latent_cov = diag(2),
    batch_means = c(-0.2, -0.1, 0),
    batch_sds = rep(0.1, 3),
    lambdas = rep(8, 3),
    thetas = rep(6, 3)
  )

  Y <- rbind(site1$Z, site2$Z)
  batch <- rep(c("a", "b"), each = 20)
  fit <- ggcombat(Y, batch, L = site1$L, max_iter = 5)

  expect_equal(dim(fit$Z_harmonized), dim(Y))
  expect_equal(dim(fit$Y_harmonized), dim(Y))
  expect_length(fit$fits, 2)
})
