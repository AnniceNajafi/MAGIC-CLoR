#' Synthetic Gaussian-mixture scRNA-seq generator with simulated dropout.
#'
#' Returns a list with:
#'   - truth: cells x genes "true" expression (positive, denoised)
#'   - obs:   cells x genes observed (sparse, with dropout)
#'   - labels: integer cluster label per cell
#'   - signal_genes / noise_genes: gene-index splits
make_synth <- function(n_cells = 300L, n_genes = 100L, n_clusters = 3L,
                       dropout = 0.6, seed = 42L) {
  set.seed(seed)
  labels <- sample.int(n_clusters, n_cells, replace = TRUE)

  # half the genes carry cluster-specific signal, half are noise.
  n_signal <- n_genes %/% 2L
  signal_genes <- seq_len(n_signal)
  noise_genes  <- setdiff(seq_len(n_genes), signal_genes)

  centroids <- matrix(0, nrow = n_clusters, ncol = n_genes)
  for (c in seq_len(n_clusters)) {
    sg <- sample(signal_genes, size = n_signal %/% n_clusters)
    centroids[c, sg] <- runif(length(sg), 2, 5)
  }
  centroids[, noise_genes] <- runif(n_clusters * length(noise_genes), 0.2, 0.6)

  truth <- centroids[labels, , drop = FALSE] +
           matrix(rnorm(n_cells * n_genes, sd = 0.15), n_cells, n_genes)
  truth[truth < 0] <- 0

  mask <- matrix(stats::rbinom(n_cells * n_genes, 1, 1 - dropout),
                 n_cells, n_genes)
  obs <- truth * mask

  list(truth = truth, obs = obs, labels = labels,
       signal_genes = signal_genes, noise_genes = noise_genes,
       n_cells = n_cells, n_genes = n_genes, n_clusters = n_clusters,
       dropout = dropout)
}

mse <- function(a, b) mean((a - b)^2)
