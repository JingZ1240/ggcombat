test_that("dense results are converted to graph-order L", {
  dense_results <- list(
    CID = c(3, 2, 1),
    Clist = c(4, 2, 1, 3, 6, 5)
  )

  graph <- ggcombat:::make_L_from_dense_results(dense_results, n_blocks = 2)

  expect_equal(dim(graph$L), c(6, 3))
  expect_equal(rowSums(graph$L), rep(1, 6))
  expect_equal(graph$CID, c(3, 2))
  expect_equal(graph$block_index, c(4, 2, 1, 3, 6))
  expect_equal(graph$singleton_index, 5)
})

test_that("Z lists can be reordered by graph Clist", {
  Z <- matrix(seq_len(12), nrow = 3)
  graph_fit <- list(results = list(Clist = c(4, 2, 1, 3)))

  out <- reorder_Z_list_by_graph(list(Z), graph_fit)

  expect_equal(out[[1]], Z[, c(4, 2, 1, 3)])
})
