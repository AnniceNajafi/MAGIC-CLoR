#' Count-based synthetic scRNA-seq for MCV tests.
#'
#' Per-cell Poisson counts around cluster centroid expression vectors. No
#' multiplicative dropout - the Poisson sampling itself creates zeros.
make_counts <- function(n_cells = 300L, n_genes = 100L, n_clusters = 3L,
                        seed = 7L, libsize = 500) {
  set.seed(seed)
  labels <- sample.int(n_clusters, n_cells, replace = TRUE)

  n_signal <- n_genes %/% 2L
  signal_genes <- seq_len(n_signal)
  noise_genes  <- setdiff(seq_len(n_genes), signal_genes)

  cent <- matrix(0.1, nrow = n_clusters, ncol = n_genes)
  for (c in seq_len(n_clusters)) {
    sg <- sample(signal_genes, size = n_signal %/% n_clusters)
    cent[c, sg] <- runif(length(sg), 3, 8)
  }
  cent[, noise_genes] <- runif(n_clusters * length(noise_genes), 0.2, 0.6)

  lam <- cent[labels, , drop = FALSE]
  lam <- lam / rowSums(lam) * libsize
  counts <- matrix(stats::rpois(n_cells * n_genes, lambda = lam),
                   nrow = n_cells, ncol = n_genes)

  list(counts = counts, lambda = lam, labels = labels,
       signal_genes = signal_genes, noise_genes = noise_genes)
}
