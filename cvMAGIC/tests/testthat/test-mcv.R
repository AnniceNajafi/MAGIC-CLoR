test_that("poisson_split sums back to original counts", {
  set.seed(11)
  X <- matrix(rpois(500, lambda = 5), 50, 10)
  sp <- cvMAGIC:::poisson_split(X, split_p = 0.5)
  expect_equal(sp$A + sp$B, X)
  expect_true(all(sp$A >= 0))
  expect_true(all(sp$B >= 0))
})

test_that("mcv_select_t errors on non-count input", {
  X <- matrix(runif(200), 20, 10)
  expect_error(mcv_select_t(X, t_values = 1:3),
               regexp = "non-negative integer count matrix")
})

test_that("MCV picks a finite optimum, not t=1 or t=T_max", {
  d <- make_counts(n_cells = 250L, n_genes = 80L, n_clusters = 3L, seed = 5L)
  res <- mcv_select_t(d$counts, t_values = 1:8,
                      npca = 12L, k = 20L, ka = 4L, seed = 5L)
  cat(sprintf("\n  MCV loss curve: %s\n",
              paste(sprintf("t=%d:%.2f", res$t_values, res$loss),
                    collapse = "  ")))
  cat(sprintf("  MCV-chosen t = %d\n", res$t))
  expect_true(res$t %in% res$t_values)
  expect_gte(res$t, 2L)
})

test_that("MCV-chosen t beats t=1 in true reconstruction (sqrt-scale)", {
  d <- make_counts(n_cells = 250L, n_genes = 80L, n_clusters = 3L, seed = 6L)
  res <- mcv_select_t(d$counts, t_values = 1:8,
                      npca = 12L, k = 20L, ka = 4L, seed = 6L)
  X_pre <- sqrt(d$counts)
  truth_pre <- sqrt(d$lambda)
  g_full <- magic_graph(X_pre, npca = 12L, k = 20L, ka = 4L)
  imp_best <- magic_impute(X_pre, t = res$t, graph = g_full)
  imp_t1   <- magic_impute(X_pre, t = 1L,   graph = g_full)
  err_best <- mean((imp_best - truth_pre)^2)
  err_t1   <- mean((imp_t1   - truth_pre)^2)
  cat(sprintf("\n  truth-MSE  t=1: %.4f   t=%d (MCV): %.4f\n",
              err_t1, res$t, err_best))
  expect_lt(err_best, err_t1)
})
