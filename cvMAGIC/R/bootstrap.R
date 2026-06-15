#' Bootstrap-based per-(cell, gene) standard error for MAGIC imputation.
#'
#' Re-runs MAGIC on `n_boot` cell bootstraps and aggregates per-(cell, gene)
#' imputed values to produce point-wise standard errors. Each bootstrap
#' resamples cells with replacement; the graph is rebuilt on the bootstrap
#' sample and imputation is mapped back to the original cell indices by
#' averaging across the duplicates a given cell received.
#'
#' This addresses the over-confidence problem documented by Andrews & Hemberg
#' (2019): MAGIC point estimates carry no native uncertainty, so downstream
#' tests assuming honest variance are mis-calibrated.
#'
#' @param X cells x genes data matrix (pre-transformed). For count data, pass
#'   `sqrt(X)`.
#' @param t Diffusion time.
#' @param n_boot Number of bootstrap replicates.
#' @param npca,k,ka Graph parameters.
#' @param seed RNG seed.
#' @return List with `mean` (cells x genes; bootstrap-mean imputed value),
#'   `se` (cells x genes; bootstrap SE), and `n_eff` (cells; number of
#'   bootstraps that included each original cell).
#' @export
magic_bootstrap_se <- function(X, t = 3L, n_boot = 30L,
                               npca = 20L, k = 30L, ka = 5L,
                               seed = 1L) {
  N <- nrow(X); G <- ncol(X)
  set.seed(seed)

  sum_x  <- matrix(0, N, G)
  sum_x2 <- matrix(0, N, G)
  n_eff  <- integer(N)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(N, N, replace = TRUE)
    Xb <- X[idx, , drop = FALSE]

    gb <- magic_graph(Xb, npca = npca, k = k, ka = ka)
    Yb <- diffuse(gb$M, Xb, t)

    for (orig in unique(idx)) {
      pos <- which(idx == orig)
      cellval <- colMeans(Yb[pos, , drop = FALSE])
      sum_x[orig, ]  <- sum_x[orig, ]  + cellval
      sum_x2[orig, ] <- sum_x2[orig, ] + cellval^2
      n_eff[orig]    <- n_eff[orig] + 1L
    }
  }

  mean_x <- sum_x / pmax(n_eff, 1L)
  var_x  <- sum_x2 / pmax(n_eff, 1L) - mean_x^2
  var_x[var_x < 0] <- 0
  se_x <- sqrt(var_x)

  unseen <- n_eff == 0L
  if (any(unseen)) {
    mean_x[unseen, ] <- NA_real_
    se_x[unseen, ]   <- NA_real_
  }

  list(mean = mean_x, se = se_x, n_eff = n_eff, n_boot = n_boot, t = t)
}
