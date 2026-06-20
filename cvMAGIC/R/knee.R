#' Knee-point selection of MAGIC's diffusion time (van Dijk et al., Fig S1C).
#'
#' Diffusion-time selection by the relative Frobenius change between
#' successive imputations,
#' delta(t) = || M^t X - M^(t-1) X ||_F / || M^(t-1) X ||_F, picking the
#' t at which delta(t) drops below a small fractional threshold (the "knee").
#' This is a rotation-free simplification of the Procrustes disparity used by
#' the upstream MAGIC auto-selection; it serves as cvMAGIC's reference
#' baseline that MCV-based selection replaces.
#'
#' @param X cells x genes data matrix (pre-transformed, e.g. sqrt(counts)).
#' @param t_max Maximum diffusion time to consider.
#' @param tol Stop when relative change drops below this value.
#' @param npca,k,ka Graph parameters.
#' @return List with `t` (chosen), `delta` (vector of relative changes), and
#'   `graph`.
#' @export
knee_select_t <- function(X, t_max = 12L, tol = 0.05,
                          npca = 20L, k = 30L, ka = 5L) {
  g <- magic_graph(X, npca = npca, k = k, ka = ka)
  prev <- as.matrix(X)
  deltas <- numeric(t_max)
  chosen <- t_max
  for (tt in seq_len(t_max)) {
    cur <- as.matrix(g$M %*% prev)
    num <- sqrt(sum((cur - prev)^2))
    den <- sqrt(sum(prev^2)) + 1e-12
    deltas[tt] <- num / den
    if (tt >= 2L && deltas[tt] < tol) {
      chosen <- tt
      break
    }
    prev <- cur
  }
  list(t = chosen, delta = deltas[seq_len(min(tt, t_max))], graph = g)
}
