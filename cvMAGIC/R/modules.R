#' MCV losses per gene across a diffusion-time sweep.
#'
#' Like `mcv_select_t` but returns a t x gene loss matrix rather than the
#' aggregate. Used internally by `magic_modules`.
#'
#' @inheritParams mcv_select_t
#' @return List with `graph`, `t_values`, `loss_per_gene` (rows = t_values,
#'   cols = genes).
#' @keywords internal
mcv_per_gene <- function(X, t_values = 1:10,
                         npca = 20L, k = 30L, ka = NULL, decay = 1,
                         split_p = 0.5, seed = 1L) {
  if (!is.numeric(X) || any(X < 0) || any(X != round(X))) {
    stop("mcv_per_gene requires a non-negative integer count matrix.")
  }
  t_values <- validate_t_values(t_values)
  set.seed(seed)
  sp <- poisson_split(X, split_p)
  A <- sp$A; B <- sp$B
  A_pre <- sqrt(A)
  g <- magic_graph(A_pre, npca = npca, k = k, ka = ka, decay = decay)

  loss_per_gene <- matrix(NA_real_, nrow = length(t_values), ncol = ncol(X))
  scale_ab <- (1 - split_p) / split_p

  Y <- A_pre
  prev_t <- 0L
  for (i in seq_along(t_values)) {
    while (prev_t < t_values[i]) {
      Y <- g$M %*% Y
      prev_t <- prev_t + 1L
    }
    Y_counts <- (as.matrix(Y))^2
    loss_per_gene[i, ] <- colMeans((Y_counts * scale_ab - B)^2)
  }

  list(graph = g, t_values = t_values, loss_per_gene = loss_per_gene,
       split_p = split_p)
}

#' MAGIC with per-gene-module adaptive diffusion time.
#'
#' Clusters genes into `n_modules` modules (k-means on gene expression
#' profiles), selects a diffusion time per module by molecular cross-validation,
#' and applies the module-specific t to its genes for the final imputation.
#'
#' The final imputation is computed directly on the input count scale (no
#' sqrt round-trip) so that the output is comparable to other MAGIC
#' implementations. The MCV split internally sqrt-transforms the Poisson
#' halves for graph construction only.
#'
#' @param X Integer count matrix, cells x genes.
#' @param n_modules Number of gene modules.
#' @param t_values Candidate diffusion times.
#' @param npca,k,ka,decay Graph parameters passed to `magic_graph`.
#' @param split_p Train-half mass for Poisson thinning.
#' @param seed RNG seed.
#' @return List with `imputed` (cells x genes, on the input count scale),
#'   `modules` (gene -> module index), `t_per_module`, `t_values`,
#'   `loss_per_gene`, and `graph`.
#' @export
magic_modules <- function(X, n_modules = 10L, t_values = 1:10,
                          npca = 20L, k = 30L, ka = NULL, decay = 1,
                          split_p = 0.5, seed = 1L) {
  n_modules <- as.integer(n_modules)
  mcv <- mcv_per_gene(X, t_values = t_values,
                      npca = npca, k = k, ka = ka, decay = decay,
                      split_p = split_p, seed = seed)

  gene_profiles <- t(scale(X, center = TRUE, scale = FALSE))
  n_pc <- min(10L, n_modules * 2L, nrow(gene_profiles) - 1L,
              ncol(gene_profiles) - 1L)
  pcs_g <- irlba::prcomp_irlba(gene_profiles, n = n_pc,
                               center = TRUE, scale. = FALSE)$x
  set.seed(seed)
  km <- stats::kmeans(pcs_g, centers = n_modules, nstart = 10L)
  modules <- km$cluster

  t_per_module <- integer(n_modules)
  for (m in seq_len(n_modules)) {
    cols <- which(modules == m)
    mloss <- rowMeans(mcv$loss_per_gene[, cols, drop = FALSE])
    t_per_module[m] <- mcv$t_values[which.min(mloss)]
  }

  g_full <- magic_graph(X, npca = npca, k = k, ka = ka, decay = decay)
  t_unique <- sort(unique(t_per_module))
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

  imputed <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X))
  for (m in seq_len(n_modules)) {
    cols <- which(modules == m); tm <- t_per_module[m]
    imputed[, cols] <- cache[[as.character(tm)]][, cols, drop = FALSE]
  }

  list(imputed = imputed,
       modules = modules,
       t_per_module = t_per_module,
       t_values = mcv$t_values,
       loss_per_gene = mcv$loss_per_gene,
       graph = g_full)
}
