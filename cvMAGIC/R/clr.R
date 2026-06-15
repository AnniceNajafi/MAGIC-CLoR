#' Context Likelihood of Relatedness (CLR) calibration of a correlation matrix.
#'
#' Raw gene-gene correlations after graph-diffusion imputation are
#' confounded: diffusion shares cell-state, lineage, and library-size
#' structure across neighbours, so many gene pairs co-vary for reasons
#' unrelated to direct regulation, and highly connected ("hub") genes
#' correlate with almost everything. The CLR transform (Faith et al.,
#' \emph{PLoS Biology} 2007) null-calibrates each association against the
#' background distribution of each gene's correlations:
#' \deqn{z_i(j) = \max\!\big(0, (C_{ij} - \mu_i)/\sigma_i\big), \quad
#'       \mathrm{CLR}_{ij} = \sqrt{z_i(j)^2 + z_j(i)^2},}
#' where \eqn{\mu_i, \sigma_i} are the mean and s.d. of gene \eqn{i}'s row
#' of \eqn{C}. A pair scores highly only if the two genes are more
#' correlated than each is \emph{typically} correlated with all genes.
#'
#' In the cvMAGIC benchmarks, replacing raw \eqn{|\mathrm{corr}|} with
#' \code{clr_from_cor()} consistently and substantially improves
#' TF-target recovery across four datasets, and the gain survives
#' detection/expression-matched negative controls. See
#' \code{\link{magic_clr}} for the end-to-end wrapper.
#'
#' @param C A symmetric gene-by-gene correlation matrix (e.g. from
#'   \code{stats::cor}). \code{NA}s are treated as 0.
#' @return A gene-by-gene CLR association matrix (non-negative, symmetric,
#'   same dimnames as \code{C}).
#' @export
clr_from_cor <- function(C) {
  C <- as.matrix(C)
  C[is.na(C)] <- 0
  mu  <- rowMeans(C)
  sdv <- apply(C, 1L, stats::sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  # Column-major recycling makes (C - mu)/sdv operate row-wise:
  # element [i, j] becomes (C[i, j] - mu[i]) / sdv[i].
  Z  <- (C - mu) / sdv
  Zp <- pmax(Z, 0)
  clr <- sqrt(Zp^2 + t(Zp)^2)
  dimnames(clr) <- dimnames(C)
  clr
}

#' MAGIC imputation with a CLR-calibrated gene-interaction readout.
#'
#' Runs MAGIC diffusion at time \code{t} (use \code{t = 0} for no
#' imputation), computes the gene-gene correlation of the imputed matrix,
#' and returns both the raw correlation and its CLR-calibrated counterpart
#' (\code{\link{clr_from_cor}}). The cvMAGIC investigation found that the
#' dominant limitation for relationship recovery is not the diffusion time
#' or graph but the \emph{raw-correlation readout}; the CLR association is
#' the recommended score for ranking gene-gene (e.g. TF-target)
#' relationships.
#'
#' @param X Integer count matrix, cells x genes (raw UMI counts).
#' @param t Diffusion time. \code{t = 0} skips imputation (CLR on raw
#'   counts), which already captures much of the gain; larger \code{t}
#'   adds value on datasets where diffusion helps.
#' @param npca,k,ka,decay Graph parameters passed to \code{\link{magic_graph}}.
#' @param method Correlation method, passed to \code{stats::cor}
#'   (default \code{"spearman"}).
#' @param graph Optional precomputed graph from \code{\link{magic_graph}}.
#' @return List with \code{imputed} (cells x genes), \code{cor} (raw
#'   gene-gene correlation), \code{assoc} (CLR-calibrated association),
#'   and \code{t}.
#' @seealso \code{\link{clr_from_cor}}
#' @export
magic_clr <- function(X, t = 3L, npca = 20L, k = 30L, ka = NULL, decay = 1,
                      method = "spearman", graph = NULL) {
  t <- as.integer(t)
  if (t == 0L) {
    imp <- as.matrix(X)
  } else {
    if (is.null(graph)) graph <- magic_graph(X, npca = npca, k = k, ka = ka,
                                             decay = decay)
    imp <- magic_impute(X, t = t, graph = graph)
  }
  colnames(imp) <- colnames(X)
  C <- suppressWarnings(stats::cor(imp, method = method))
  C[is.na(C)] <- 0
  list(imputed = imp, cor = C, assoc = clr_from_cor(C), t = t)
}
