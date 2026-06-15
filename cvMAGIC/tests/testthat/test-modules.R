test_that("magic_modules picks a t per module and applies it", {
  d <- make_mixed_counts(n_cells = 250L, n_genes = 90L,
                         n_clusters = 3L, n_groups = 3L,
                         scales = c(0.3, 1.5, 6.0), seed = 21L)
  out <- magic_modules(d$counts, n_modules = 4L, t_values = 1:8,
                       npca = 12L, k = 20L, seed = 21L)
  expect_equal(dim(out$imputed), dim(d$counts))
  expect_length(out$t_per_module, 4L)
  expect_true(all(out$t_per_module %in% out$t_values))
  cat(sprintf("\n  t_per_module = %s\n",
              paste(out$t_per_module, collapse = ", ")))
})

test_that("per-module t is no worse than the single best global t on mixed data", {
  d <- make_mixed_counts(n_cells = 300L, n_genes = 120L,
                         n_clusters = 3L, n_groups = 3L,
                         scales = c(0.3, 1.5, 6.0), seed = 22L)
  truth <- d$lambda

  out <- magic_modules(d$counts, n_modules = 6L, t_values = 1:8,
                       npca = 14L, k = 22L, seed = 22L)
  err_mod <- mean((out$imputed - truth)^2)

  X_in <- d$counts
  g_full <- out$graph
  errs_global <- numeric(length(out$t_values))
  Y <- X_in; prev_t <- 0L
  for (i in seq_along(out$t_values)) {
    while (prev_t < out$t_values[i]) {
      Y <- g_full$M %*% Y; prev_t <- prev_t + 1L
    }
    errs_global[i] <- mean((as.matrix(Y) - truth)^2)
  }
  best_global_err <- min(errs_global)
  best_global_t   <- out$t_values[which.min(errs_global)]

  cat(sprintf("\n  best global t = %d, err = %.5f\n",
              best_global_t, best_global_err))
  cat(sprintf("  per-module err = %.5f (t_per_module = %s)\n",
              err_mod, paste(out$t_per_module, collapse = ",")))
  expect_lte(err_mod, best_global_err * 1.02)
})
