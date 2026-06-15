#' Build the MAGIC cell-cell Markov affinity matrix.
#'
#' Performs PCA preprocessing, finds k nearest neighbours in PC space, applies
#' an adaptive alpha-decay kernel (default exponential, matching the upstream
#' Python MAGIC implementation), symmetrises, and row-normalises to produce a
#' sparse Markov transition matrix.
#'
#' Defaults track the upstream MAGIC implementation: decay = 1 (exponential
#' kernel) and bandwidth = distance to the k-th (i.e. last) neighbour. The
#' published paper describes a Gaussian kernel with bandwidth set to the
#' ka-th neighbour; pass `decay = 2` and `ka < k` to use that variant.
#'
#' @param X Numeric matrix of cells (rows) by genes (columns). Should be
#'   library-size normalised and sqrt or log transformed before calling.
#' @param npca Number of principal components to use for the kNN search.
#' @param k Number of nearest neighbours per cell.
#' @param ka Bandwidth neighbour index. If `NULL` (default) use `k`, matching
#'   the upstream implementation. Otherwise the kernel width is the distance
#'   to the ka-th neighbour (paper's spec).
#' @param decay Alpha-decay kernel exponent: `K(i,j) = exp(-(d/sigma)^decay)`.
#'   Default 1 (exponential) matches upstream; set 2 for Gaussian.
#' @return A list with `M` (sparse N x N Markov matrix), `pca` (cells x npca),
#'   and the parameters used.
#' @export
magic_graph <- function(X, npca = 20L, k = 30L, ka = NULL, decay = 1) {
  if (!is.matrix(X) && !inherits(X, "Matrix")) {
    stop("X must be a matrix (cells x genes).")
  }
  N <- nrow(X)
  if (k >= N) stop("k must be < number of cells.")
  if (is.null(ka)) ka <- k
  if (ka > k)  stop("ka must be <= k.")
  if (decay <= 0) stop("decay must be > 0.")

  npca <- min(npca, ncol(X) - 1L, N - 1L)
  pcs <- irlba::prcomp_irlba(X, n = npca, center = TRUE, scale. = FALSE)$x

  nn <- RANN::nn2(pcs, pcs, k = k + 1L)
  idx <- nn$nn.idx[, -1, drop = FALSE]
  dst <- nn$nn.dists[, -1, drop = FALSE]

  sigma <- dst[, ka]
  sigma[sigma == 0] <- .Machine$double.eps

  aff <- exp(-(dst / sigma)^decay)

  rows <- rep.int(seq_len(N), times = k)
  cols <- as.integer(idx)
  vals <- as.numeric(aff)

  A <- Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(N, N))
  A <- (A + Matrix::t(A)) / 2

  rs <- Matrix::rowSums(A)
  rs[rs == 0] <- 1
  M <- A / rs

  list(M = M, pca = pcs, npca = npca, k = k, ka = ka, decay = decay)
}

#' Apply diffusion: compute the action of M^t on X.
#'
#' Computes M to the power t multiplied by X without forming a dense M^t,
#' by iterating matrix-vector multiplications.
#'
#' @param M Sparse Markov matrix (cells x cells).
#' @param X Original count / normalised data matrix (cells x genes).
#' @param t Diffusion time (non-negative integer).
#' @return Imputed matrix of the same shape as X.
#' @keywords internal
diffuse <- function(M, X, t) {
  if (t < 0) stop("t must be >= 0.")
  if (t == 0) return(as.matrix(X))
  Y <- X
  for (i in seq_len(t)) Y <- M %*% Y
  as.matrix(Y)
}
