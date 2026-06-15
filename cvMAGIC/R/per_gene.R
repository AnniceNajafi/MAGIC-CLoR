#' MAGIC with per-gene adaptive diffusion time.
#'
#' For each gene, select the diffusion time that minimises its
#' per-gene molecular cross-validation loss, then assemble the imputed
#' matrix column-wise using each gene's chosen \eqn{t^\star_g}.
#'
#' This is the limit case of `magic_modules` with `n_modules = ncol(X)`:
#' the k-means step is removed and each gene is its own module. Where
#' `magic_modules` averages per-gene MCV losses within a cluster before
#' picking a single t, this function preserves per-gene granularity.
#'
#' Two MCV loss flavours are supported:
#'   * `"mse"` (default): the original count-scale MSE used by
#'     [mcv_select_t()].
#'   * `"pearson"`: per-(cell, gene) Pearson residuals - the squared
#'     standardised residual `(B - prediction)^2 / max(prediction, eps)`.
#'     This is scale-invariant per gene and gives the optimisation a
#'     less heteroscedastic loss surface.
#'
#' @param X Integer count matrix, cells x genes.
#' @param t_values Candidate diffusion times.
#' @param loss Either `"mse"` or `"pearson"`.
#' @param npca,k,ka,decay Graph parameters passed to `magic_graph`.
#' @param split_p Train-half mass for Poisson thinning.
#' @param seed RNG seed.
#' @param eps Floor for the Pearson denominator (avoids division by zero
#'   on gene-cell pairs that the train graph predicts as zero).
#' @return List with `imputed` (cells x genes, count scale),
#'   `t_per_gene` (length-G integer), `t_values`,
#'   `loss_per_gene` (length(t_values) x G), `graph`, and the chosen
#'   `loss` flavour.
#' @export
magic_per_gene <- function(X, t_values = 1:10,
                           loss = c("mse", "pearson"),
                           npca = 20L, k = 30L, ka = NULL, decay = 1,
                           split_p = 0.5, seed = 1L,
                           eps = 1e-6) {
  loss <- match.arg(loss)
  if (!is.numeric(X) || any(X < 0) || any(X != round(X))) {
    stop("magic_per_gene requires a non-negative integer count matrix.")
  }
  set.seed(seed)
  sp <- poisson_split(X, split_p)
  A <- sp$A; B <- sp$B
  A_pre <- sqrt(A)
  g_train <- magic_graph(A_pre, npca = npca, k = k, ka = ka, decay = decay)

  t_values <- sort(unique(as.integer(t_values)))
  G <- ncol(X)
  loss_per_gene <- matrix(NA_real_, nrow = length(t_values), ncol = G)
  scale_ab <- (1 - split_p) / split_p

  Y <- A_pre
  prev_t <- 0L
  for (i in seq_along(t_values)) {
    while (prev_t < t_values[i]) {
      Y <- g_train$M %*% Y
      prev_t <- prev_t + 1L
    }
    Y_counts <- (as.matrix(Y))^2
    predicted <- Y_counts * scale_ab
    if (loss == "mse") {
      loss_per_gene[i, ] <- colMeans((predicted - B)^2)
    } else {
      denom <- pmax(predicted, eps)
      loss_per_gene[i, ] <- colMeans((B - predicted)^2 / denom)
    }
  }

  t_per_gene <- t_values[apply(loss_per_gene, 2, which.min)]

  g_full <- magic_graph(X, npca = npca, k = k, ka = ka, decay = decay)
  t_unique <- sort(unique(t_per_gene))
  cache <- vector("list", length(t_unique))
  names(cache) <- as.character(t_unique)
  Y <- X
  prev_t <- 0L
  for (tu in t_unique) {
    while (prev_t < tu) {
      Y <- g_full$M %*% Y
      prev_t <- prev_t + 1L
    }
    cache[[as.character(tu)]] <- as.matrix(Y)
  }

  imputed <- matrix(NA_real_, nrow = nrow(X), ncol = G)
  for (g in seq_len(G)) {
    imputed[, g] <- cache[[as.character(t_per_gene[g])]][, g]
  }

  list(imputed = imputed,
       t_per_gene = t_per_gene,
       t_values = t_values,
       loss_per_gene = loss_per_gene,
       graph = g_full,
       loss = loss)
}
