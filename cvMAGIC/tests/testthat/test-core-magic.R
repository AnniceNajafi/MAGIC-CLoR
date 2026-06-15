test_that("magic_graph returns a row-stochastic sparse matrix", {
  d <- make_synth(n_cells = 150L, n_genes = 60L, n_clusters = 3L,
                  dropout = 0.5, seed = 1L)
  g <- magic_graph(d$obs, npca = 10L, k = 15L, ka = 4L)
  expect_s4_class(g$M, "Matrix")
  expect_equal(nrow(g$M), nrow(d$obs))
  rs <- Matrix::rowSums(g$M)
  expect_true(all(abs(rs - 1) < 1e-8))
})

test_that("MAGIC reduces MSE vs raw on synthetic dropout data", {
  d <- make_synth(n_cells = 250L, n_genes = 80L, n_clusters = 3L,
                  dropout = 0.6, seed = 2L)
  imp <- magic(d$obs, t = 3L, npca = 15L, k = 20L, ka = 4L)
  err_raw   <- mse(d$obs,  d$truth)
  err_magic <- mse(imp,    d$truth)
  cat(sprintf("\n  raw MSE = %.4f, MAGIC MSE = %.4f\n", err_raw, err_magic))
  expect_lt(err_magic, err_raw)
})

test_that("Imputation preserves cell count and gene count", {
  d <- make_synth(n_cells = 120L, n_genes = 40L, dropout = 0.5, seed = 3L)
  imp <- magic(d$obs, t = 2L, npca = 8L, k = 12L, ka = 3L)
  expect_equal(dim(imp), dim(d$obs))
})

test_that("diffuse(t=0) is identity on X", {
  d <- make_synth(n_cells = 80L, n_genes = 30L, seed = 4L)
  g <- magic_graph(d$obs, npca = 6L, k = 10L, ka = 3L)
  Y <- cvMAGIC:::diffuse(g$M, d$obs, t = 0L)
  expect_equal(as.matrix(Y), as.matrix(d$obs))
})
