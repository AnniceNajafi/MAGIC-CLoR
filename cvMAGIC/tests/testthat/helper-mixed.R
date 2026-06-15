#' Mixed-scale synthetic counts: gene groups with very different library shares.
#'
#' Creates `n_groups` equal-sized gene groups, each with a different mean
#' expression scale. This produces modules where the optimal MAGIC t differs:
#' low-count groups need more diffusion to overcome zeros; high-count groups
#' need less.
make_mixed_counts <- function(n_cells = 300L, n_genes = 120L,
                              n_clusters = 3L, n_groups = 3L,
                              scales = c(0.3, 1.5, 6.0),
                              seed = 13L, libsize = 1500) {
  set.seed(seed)
  stopifnot(length(scales) == n_groups, n_genes %% n_groups == 0)
  per_group <- n_genes %/% n_groups
  group <- rep(seq_len(n_groups), each = per_group)

  labels <- sample.int(n_clusters, n_cells, replace = TRUE)
  cent <- matrix(0.1, nrow = n_clusters, ncol = n_genes)
  for (g in seq_len(n_groups)) {
    gcols <- which(group == g)
    for (c in seq_len(n_clusters)) {
      sg <- sample(gcols, size = per_group %/% n_clusters)
      cent[c, sg] <- runif(length(sg), 2, 5) * scales[g]
    }
  }

  lam <- cent[labels, , drop = FALSE]
  lam <- lam / rowSums(lam) * libsize
  counts <- matrix(stats::rpois(n_cells * n_genes, lambda = lam),
                   nrow = n_cells, ncol = n_genes)

  list(counts = counts, lambda = lam, labels = labels,
       group = group, n_groups = n_groups, scales = scales)
}
