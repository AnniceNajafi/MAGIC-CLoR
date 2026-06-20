#' Run MAGIC imputation with a fixed diffusion time.
#'
#' @param X Cells x genes data matrix (library-size normalised + sqrt-transformed
#'   is recommended). Must be the same matrix used for graph construction OR a
#'   compatible matrix of identical row count.
#' @param t Diffusion time (non-negative integer; default 3).
#' @param graph Optional pre-computed graph from `magic_graph()`. If NULL, one is
#'   built from X.
#' @param ... Passed to `magic_graph()` if `graph` is NULL.
#' @return Imputed cells x genes matrix.
#' @export
magic_impute <- function(X, t = 3L, graph = NULL, ...) {
  if (is.null(graph)) graph <- magic_graph(X, ...)
  if (nrow(graph$M) != nrow(X)) {
    stop("graph$M and X must have the same number of cells.")
  }
  diffuse(graph$M, X, t)
}

#' Build the graph and run MAGIC in one call, with a fixed t.
#' @inheritParams magic_impute
#' @param npca Number of principal components for the kNN search.
#' @param k Number of nearest neighbours per cell.
#' @param ka Adaptive bandwidth parameter (kernel width = distance to ka-th
#'   neighbour).
#' @return Imputed cells x genes matrix.
#' @export
magic <- function(X, t = 3L, npca = 20L, k = 30L, ka = 5L) {
  g <- magic_graph(X, npca = npca, k = k, ka = ka)
  magic_impute(X, t = t, graph = g)
}
