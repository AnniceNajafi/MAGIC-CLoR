#' Validate a Poisson-thinning split probability.
#'
#' `split_p` must be a single finite number strictly inside (0, 1): at
#' `split_p = 0` the rescale `(1 - split_p) / split_p` divides by zero, and
#' at `split_p = 1` the held-out half `B` is all zeros so the MCV loss
#' collapses independently of prediction quality. Either boundary (or an
#' out-of-range value) would let a selector report a confident `t` from an
#' invalid cross-validation experiment instead of failing fast.
#'
#' @param split_p Candidate train-half mass.
#' @return `split_p` invisibly, if valid; otherwise an error.
#' @keywords internal
validate_split_p <- function(split_p) {
  if (!is.numeric(split_p) || length(split_p) != 1L ||
      !is.finite(split_p) || split_p <= 0 || split_p >= 1) {
    stop("split_p must be a single finite number strictly between 0 and 1 ",
         "(got ", paste(format(split_p), collapse = ", "), ").")
  }
  invisible(split_p)
}

#' Validate and canonicalise a diffusion-time search grid.
#'
#' Candidate diffusion times must be non-negative integers (`t = 0` is the
#' valid "no imputation" option). Negative, non-integer, `NA`, or empty
#' grids are rejected rather than silently coerced, since the forward sweep
#' assumes candidates are ordered up from `prev_t = 0`: a negative candidate
#' would skip the loop entirely and leave the loss attached to an impossible
#' `t`. Returns the sorted, de-duplicated integer grid.
#'
#' @param t_values Candidate diffusion times.
#' @return Sorted unique integer vector of valid diffusion times.
#' @keywords internal
validate_t_values <- function(t_values) {
  if (!is.numeric(t_values) || length(t_values) == 0L ||
      any(!is.finite(t_values))) {
    stop("t_values must be a non-empty numeric vector of finite, ",
         "non-negative integers.")
  }
  if (any(t_values < 0) || any(t_values != round(t_values))) {
    stop("t_values must be non-negative integers (got ",
         paste(format(t_values), collapse = ", "), ").")
  }
  sort(unique(as.integer(round(t_values))))
}

#' Poisson-thinning split for molecular cross-validation (Batson et al. 2019).
#'
#' For an integer count matrix X, returns two matrices A and B with
#' A_ij ~ Binomial(X_ij, split_p) and B_ij = X_ij - A_ij. Under the Poisson
#' assumption, A and B are independent with means split_p * lambda and
#' (1 - split_p) * lambda respectively.
#'
#' @param X Integer (count) matrix, cells x genes.
#' @param split_p Probability mass for the training half.
#' @return List with `A` and `B` matrices.
#' @keywords internal
poisson_split <- function(X, split_p = 0.5) {
  validate_split_p(split_p)
  X <- as.matrix(X)
  storage.mode(X) <- "integer"
  A <- matrix(stats::rbinom(length(X), X, split_p),
              nrow = nrow(X), ncol = ncol(X))
  B <- X - A
  list(A = A, B = B)
}

#' Select MAGIC's diffusion time via molecular cross-validation.
#'
#' Splits a count matrix into two independent Poisson halves, builds the
#' MAGIC graph on the training half, sweeps diffusion times, and picks
#' the t that minimises validation MSE against the held-out half (with
#' the correct scaling so both halves estimate comparable quantities).
#'
#' @param X Integer count matrix, cells x genes.
#' @param t_values Candidate diffusion times (non-negative integers).
#'   Default `0:10` - `t = 0` represents "no imputation", letting MCV
#'   reject smoothing on datasets where the raw counts already carry
#'   the structure of interest (e.g. low-dropout 10x), letting MCV reject
#'   smoothing where the raw counts already carry the signal.
#' @param npca,k,ka Graph parameters passed to `magic_graph`.
#' @param split_p Train-half mass for Poisson thinning.
#' @param seed RNG seed for reproducibility.
#' @return List with `t` (chosen), `t_values`, `loss` (per t),
#'   and `graph` (the graph built on the training half).
#' @export
mcv_select_t <- function(X, t_values = 0:10,
                         npca = 20L, k = 30L, ka = 5L,
                         split_p = 0.5, seed = 1L) {
  if (!is.numeric(X) || any(X < 0) || any(X != round(X))) {
    stop("mcv_select_t requires a non-negative integer count matrix. ",
         "Apply MCV before normalisation/sqrt transform.")
  }
  t_values <- validate_t_values(t_values)
  set.seed(seed)
  sp <- poisson_split(X, split_p = split_p)
  A <- sp$A; B <- sp$B

  A_pre <- sqrt(A)
  g <- magic_graph(A_pre, npca = npca, k = k, ka = ka)

  losses <- numeric(length(t_values))
  scale_ab <- (1 - split_p) / split_p

  Y <- A_pre
  prev_t <- 0L
  for (i in seq_along(t_values)) {
    while (prev_t < t_values[i]) {
      Y <- g$M %*% Y
      prev_t <- prev_t + 1L
    }
    Y_counts <- (as.matrix(Y))^2
    losses[i] <- mean((Y_counts * scale_ab - B)^2)
  }

  list(t = t_values[which.min(losses)],
       t_values = t_values,
       loss = losses,
       graph = g,
       split_p = split_p)
}

#' Select MAGIC's diffusion time via tolerance-regularised MCV.
#'
#' Standard `mcv_select_t` picks the argmin of the Poisson cross-validation
#' loss curve. Empirically that under-smooths sparse data (e.g. MARS-seq)
#' for the purpose of gene-gene relationship recovery: the AUC of TF-target
#' predictions keeps improving for many *t* beyond MCV's argmin, while the
#' MCV loss grows only slowly. Conversely on dense data (e.g. 10x), the
#' argmin is at *t = 0* (no smoothing) and any increase in t hurts both
#' MCV and AUC sharply.
#'
#' This selector picks the **largest** *t* whose MCV loss is within a
#' fractional tolerance `tolerance` of the minimum:
#'
#' `t* = max { t : MCV(t) <= MCV(t_argmin) * (1 + tolerance) }`
#'
#' On sparse trajectory data the rule can push to slightly larger t than the
#' argmin, but in our benchmarks the resulting AUC change is within sampling
#' noise and does not generalise across datasets (see the accompanying paper);
#' it is kept exported as a documented negative result. On dense 10x data the
#' loss curve rises steeply from t=0, so the selector stays at t=0. The default
#' `tolerance = 0.05` does not regress on dense data.
#'
#' Conceptually this is the inverse-direction sibling of the classical
#' "1-SE rule" in CV (Hastie, Tibshirani & Friedman 2009, section 7.10): the
#' 1-SE rule picks the *smallest* model whose CV loss is within 1 SE of
#' the best, on the premise that simpler is better. For graph-diffusion
#' imputation, **larger** *t* is the more regularised (more averaging
#' per cell) choice for gene-gene relationship recovery, and the
#' tolerance is fractional rather than SE-scaled because the empirical
#' MCV-loss SE is tiny relative to the AUC-relevant t range.
#'
#' @param X Integer count matrix, cells x genes.
#' @param t_values Candidate diffusion times (default `0:10`).
#' @param tolerance Fractional tolerance above the MCV minimum (default
#'   `0.05`). Set to `0` to recover the vanilla `mcv_select_t` argmin.
#' @param npca,k,ka Graph parameters passed to `magic_graph`.
#' @param split_p Train-half mass for Poisson thinning.
#' @param seed RNG seed for reproducibility.
#' @return List with `t` (chosen via tolerance rule), `t_argmin` (vanilla
#'   MCV's choice for comparison), `t_values`, `loss`, `threshold`,
#'   `tolerance`, `graph`, `split_p`.
#' @export
mcv_select_t_tolerant <- function(X, t_values = 0:10, tolerance = 0.05,
                                  npca = 20L, k = 30L, ka = 5L,
                                  split_p = 0.5, seed = 1L) {
  if (!is.numeric(tolerance) || tolerance < 0) {
    stop("tolerance must be a non-negative number.")
  }
  res <- mcv_select_t(X, t_values = t_values, npca = npca, k = k, ka = ka,
                      split_p = split_p, seed = seed)
  threshold <- min(res$loss) * (1 + tolerance)
  acceptable <- which(res$loss <= threshold)
  t_choice <- res$t_values[max(acceptable)]
  list(t = t_choice,
       t_argmin = res$t,
       t_values = res$t_values,
       loss = res$loss,
       threshold = threshold,
       tolerance = tolerance,
       graph = res$graph,
       split_p = split_p)
}

#' MAGIC with tolerance-regularised MCV t-selection.
#'
#' Picks `t` with [mcv_select_t_tolerant()], then runs [magic_impute()] with
#' that t on a graph rebuilt from the full data. This is the
#' tolerant MAGIC workflow we use in the tolerance benchmark.
#'
#' @inheritParams mcv_select_t_tolerant
#' @return List with `imputed` (cells x genes), `t`, `t_argmin`,
#'   `tolerance`, `loss`, `graph`.
#' @export
magic_tolerant <- function(X, t_values = 0:10, tolerance = 0.05,
                           npca = 20L, k = 30L, ka = 5L,
                           split_p = 0.5, seed = 1L) {
  sel <- mcv_select_t_tolerant(X, t_values = t_values,
                               tolerance = tolerance, npca = npca,
                               k = k, ka = ka, split_p = split_p,
                               seed = seed)
  if (sel$t == 0L) {
    return(list(imputed = as.matrix(X), t = 0L, t_argmin = sel$t_argmin,
                tolerance = tolerance, loss = sel$loss, graph = NULL))
  }
  g_full <- magic_graph(X, npca = npca, k = k, ka = ka)
  imp <- magic_impute(X, t = sel$t, graph = g_full)
  list(imputed = imp, t = sel$t, t_argmin = sel$t_argmin,
       tolerance = tolerance, loss = sel$loss, graph = g_full)
}
