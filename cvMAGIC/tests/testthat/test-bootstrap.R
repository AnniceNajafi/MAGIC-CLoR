test_that("magic_bootstrap_se returns sensible mean and se matrices", {
  d <- make_synth(n_cells = 120L, n_genes = 40L, n_clusters = 3L,
                  dropout = 0.5, seed = 31L)
  out <- magic_bootstrap_se(d$obs, t = 3L, n_boot = 8L,
                            npca = 8L, k = 12L, ka = 3L, seed = 31L)
  expect_equal(dim(out$mean), dim(d$obs))
  expect_equal(dim(out$se),   dim(d$obs))
  expect_true(all(is.finite(out$se[out$n_eff > 0, ])))
  expect_true(all(out$se[out$n_eff > 0, ] >= 0))
  cat(sprintf("\n  bootstrap SE: median %.4f, max %.4f\n",
              median(out$se[out$n_eff > 0, ]),
              max(out$se[out$n_eff > 0, ])))
})

test_that("Bootstrap CV (SE/|mean|) is larger for sparser genes", {
  set.seed(101)
  N <- 200L
  truth <- matrix(stats::rnorm(N * 30, mean = 3, sd = 1), N, 30)
  truth[truth < 0] <- 0
  dropout_per_gene <- c(rep(0.2, 15), rep(0.85, 15))
  obs <- truth
  for (j in seq_len(30L)) {
    keep <- stats::rbinom(N, 1, 1 - dropout_per_gene[j])
    obs[, j] <- obs[, j] * keep
  }
  out <- magic_bootstrap_se(obs, t = 3L, n_boot = 10L,
                            npca = 6L, k = 12L, ka = 3L, seed = 32L)

  eps <- 1e-6
  cv <- out$se / (abs(out$mean) + eps)
  cv_dense  <- mean(cv[, 1:15], na.rm = TRUE)
  cv_sparse <- mean(cv[, 16:30], na.rm = TRUE)
  cat(sprintf("\n  mean SE: dense=%.3f sparse=%.3f\n",
              mean(out$se[, 1:15],  na.rm = TRUE),
              mean(out$se[, 16:30], na.rm = TRUE)))
  cat(sprintf("  mean CV: dense=%.3f sparse=%.3f\n", cv_dense, cv_sparse))
  expect_gt(cv_sparse, cv_dense)
})
