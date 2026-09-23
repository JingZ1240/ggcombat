test_that("unified harmonizer supports average and reference targets", {
  set.seed(2)
  s1 <- simulate_ggcombat_data(30, c(5, 5, 2), diag(2), c(0.2, 0.1, 0), rep(0.1, 3), rep(8, 3), rep(6, 3))
  s2 <- simulate_ggcombat_data(30, c(5, 5, 2), diag(2), c(-0.2, 0.2, 0), rep(0.1, 3), rep(8, 3), rep(6, 3))
  Z_list <- list(s1$Z, s2$Z)
  fits <- fit_ggcombat_params(Z_list, s1$L, max_iter = 5)

  avg <- harmonize_ggcombat_sites(Z_list, fits = fits, L = s1$L, target = "average")
  ref <- harmonize_ggcombat_sites(
    Z_list,
    fits = fits,
    L = s1$L,
    target = "reference",
    reference = 1
  )

  expect_equal(length(avg), 2)
  expect_equal(length(ref), 2)
  expect_equal(dim(avg[[1]]), dim(s1$Z))
  expect_equal(dim(ref[[2]]), dim(s2$Z))
})
