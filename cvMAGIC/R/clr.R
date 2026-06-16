#' Apply the CLR (Context Likelihood of Relatedness) transform to a correlation matrix
#'
#' After MAGIC diffusion the gene-gene correlations get inflated, the diffusion
#' mixes cell state, lineage and library size across neighbor cells so plenty of
#' pairs end up correlated without being directly related and the hub genes
#' correlate with almost everything. CLR (Faith et al, PLoS Biology 2007)
#' rescores each pair against the spread of correlations each gene has with all
#' the other genes, so a pair only scores high when the two genes are more
#' correlated with each other than either one usually is with everything else.
#'
#' @param C a symmetric gene by gene correlation matrix, NAs are set to 0
#' @return a gene by gene CLR association matrix, non-negative and symmetric
#' @export
clr_from_cor <- function(C) {
  C <- as.matrix(C)
  C[is.na(C)] <- 0
  mu  <- rowMeans(C)
  sdv <- apply(C, 1L, stats::sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  #column-major recycling makes (C-mu)/sdv go row-wise, so [i,j] becomes (C[i,j]-mu[i])/sdv[i]
  Z  <- (C - mu) / sdv
  Zp <- pmax(Z, 0)
  clr <- sqrt(Zp^2 + t(Zp)^2)
  dimnames(clr) <- dimnames(C)
  clr
}

#' Run MAGIC and score the gene pairs with the CLR readout
#'
#' Diffuses the data for t steps (set t = 0 to skip the imputation and just
#' score the raw counts), takes the gene-gene correlation of the result and
#' returns the raw correlation together with its CLR version. In our tests the
#' CLR score picks out TF-target pairs better than the plain correlation, so
#' that is the score we use for ranking gene pairs.
#'
#' @param X integer count matrix, cells x genes (raw UMI counts)
#' @param t diffusion time, t = 0 skips imputation and scores the raw counts
#' @param npca,k,ka,decay graph parameters passed to magic_graph
#' @param method correlation method passed to stats::cor (default "spearman")
#' @param graph optional precomputed graph from magic_graph
#' @return a list with imputed (cells x genes), cor (raw correlation), assoc
#'   (the CLR association) and t
#' @seealso clr_from_cor
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
