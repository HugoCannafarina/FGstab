# =============================================================================
# Stability Selection for the Fine-Gray Competing Risks Model
# Method: Integrated Path Stability Selection (IPSS)
# =============================================================================


# =============================================================================
# S3 Class: ipss_result
# =============================================================================

#' Constructor for ipss_result objects
#'
#' @param selected_features Integer vector of selected feature indices.
#' @param stability_paths Matrix of stability paths (n_alphas x p).
#' @param efp_scores Named list of EFP scores (one per feature).
#' @param q_values Named numeric vector of q-values derived from EFP scores.
#'
#' @return An object of class \code{ipss_result}.
#'
#' @keywords internal
new_ipss_result <- function(selected_features,
                            stability_paths,
                            efp_scores,
                            q_values) {
  structure(
    list(
      selected_features = selected_features,
      stability_paths   = stability_paths,
      efp_scores        = efp_scores,
      q_values          = q_values
    ),
    class = "ipss_result"
  )
}


#' Print method for ipss_result
#'
#' @param x An \code{ipss_result} object.
#' @param ... Additional arguments (ignored).
#' @importFrom fastcmprsk Crisk fastCrrp
#' @importFrom parallel mclapply detectCores
#' @export
print.ipss_result <- function(x, ...) {
  cat("=== FGstab Result (IPSS) ===\n")
  cat("Features tested  :", ncol(x$stability_paths), "\n")
  cat("Features selected:", length(x$selected_features), "\n")
  if (length(x$selected_features) > 0) {
    cat("Selected indices :", paste(x$selected_features, collapse = ", "), "\n")
  } else {
    cat("No features selected",
        "(set target_fp or target_fdr to select features).\n")
  }
  invisible(x)
}


#' Summary method for ipss_result
#'
#' @param object An \code{ipss_result} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.ipss_result <- function(object, ...) {
  cat("=== Summary: FGstab Result (IPSS) ===\n")
  cat("Features tested        :", ncol(object$stability_paths), "\n")
  cat("Regularisation steps   :", nrow(object$stability_paths), "\n")
  cat("Features selected      :", length(object$selected_features), "\n")

  cat("\nEFP scores (top 10 most stable features):\n")
  efp_vec <- sort(unlist(object$efp_scores))[1:min(10, length(object$efp_scores))]
  print(round(efp_vec, 4))

  if (length(object$selected_features) > 0) {
    cat("\nQ-values (selected features):\n")
    sel <- as.character(object$selected_features)
    print(round(object$q_values[sel], 4))
  }

  invisible(object)
}


# =============================================================================
# Main function: FGstab()
# =============================================================================

#' Variable Selection via IPSS for the Penalised Fine-Gray Model
#'
#' Implements Integrated Path Stability Selection (IPSS) (Melikechi & Miller, 2026)
#' for the Fine-Gray competing-risks regression model. The method repeatedly
#' subsamples the data, fits a LASSO-penalised Fine-Gray model on each subsample
#' across a regularisation path, and aggregates selection frequencies into
#' stability paths. EFP scores and q-values are derived from the integrated paths
#' to control either the expected number of false positives (\code{target_fp})
#' or the false discovery rate (\code{target_fdr}). When \code{y[, 2]} contains
#' only 0 and 1 (no competing event), the Fine-Gray model reduces to the Cox
#' model; \code{FGstab} then performs IPSS for penalised Cox regression.
#' See \code{\link{sim_survival}} to simulate data in that setting.
#'
#' @param X Numeric matrix of covariates, dimensions n x p. Used as-is
#'   (no internal standardisation).
#' @param y Two-column numeric matrix. The first column must contain the event
#'   time (observed follow-up time). The second column must contain the event
#'   type indicator: 0 = censored, 1 = event of interest, 2 = competing event
#'   (omit 2 for a standard survival / Cox setting).
#' @param n_alphas Integer. Number of regularisation values (lambda) on the
#'   path. Default \code{15}.
#' @param target_fp Numeric or \code{NULL}. Maximum expected number of false
#'   positives. Features with EFP score \eqn{\leq} \code{target_fp} are
#'   selected. Ignored if \code{target_fdr} is also provided. Default
#'   \code{NULL}.
#' @param target_fdr Numeric in \eqn{(0, 1)} or \code{NULL}. Maximum false
#'   discovery rate. Features with q-value \eqn{\leq} \code{target_fdr} are
#'   selected. Takes priority over \code{target_fp}. Default 0.1.
#' @param ipss_function Character string. IPSS bounding function controlling
#'   conservatism. One of \code{"h1"}, \code{"h2"} (default), or \code{"h3"}.
#'   \code{"h3"} is the least conservative.
#' @param B Integer. Number of subsampling iterations. Default \code{50}.
#' @param n_jobs Integer. Number of cores for parallel subsampling. Default
#'   \code{1}. Increase for faster computation on multi-core machines.
#'
#' @return An object of class \code{ipss_result} with fields:
#'   \describe{
#'     \item{\code{selected_features}}{Integer vector of selected feature
#'       indices. Empty if neither \code{target_fp} nor \code{target_fdr}
#'       is specified.}
#'     \item{\code{stability_paths}}{Matrix (n_alphas x p) of selection
#'       frequencies along the regularisation path.}
#'     \item{\code{efp_scores}}{Named list of EFP scores. Lower values
#'       indicate more reliably selected features.}
#'     \item{\code{q_values}}{Named numeric vector of q-values derived from
#'       EFP scores. Interpretable as an upper bound on the FDR.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' dat <- sim_competing(n_ind = 250, n_pred = 10, n_var = 500,
#'                      design = "independent")
#' y <- as.matrix(dat$surv[, c("Tobs", "evtype")])
#'
#' # Return scores only (inspect before thresholding)
#' res <- FGstab(dat$X, y, n_alphas = 15)
#' print(res)
#' plot(res)
#'
#' # Control expected false positives
#' res_fp <- FGstab(dat$X, y, n_alphas = 15, target_fp = 1)
#' print(res_fp)
#'
#' # Control FDR
#' res_fdr <- FGstab(dat$X, y, n_alphas = 15, target_fdr = 0.2)
#' print(res_fdr)
#' }
#'
#'
#' @seealso \code{\link{print.ipss_result}}, \code{\link{summary.ipss_result}},
#'   \code{\link{plot.ipss_result}}, \code{\link{sim_competing}},
#'   \code{\link{sim_survival}}
#'   
#' @references
#' \itemize{
#'   \item Melikechi, O., & Miller, J. W. (2026). Integrated path stability
#'         selection. \emph{Journal of the American Statistical Association},
#'         121(553), 454–464.
#'   \item Fine, J. P. & Gray, R. J. (1999). A proportional hazards model for
#'         the subdistribution of a competing risk. \emph{JASA}, 94(446),
#'         496–509.
#' }
#'
#' @export
FGstab <- function(X,
                   y,
                   n_alphas      = 15L,
                   B             = 50L,
                   target_fp     = NULL,
                   target_fdr    = 0.1,
                   ipss_function = "h2",
                   n_jobs        = 1L) {

  # --- Internal constants --------------------------------------------------
  delta  <- 1     # Integration weight exponent (1 = uniform on log-lambda)
  cutoff <- 0.05  # IPSS integral early-stop threshold

  p <- ncol(X)

  # --- Regularisation path -------------------------------------------------
  alphas   <- .compute_alphas(X, y, n_alphas, rep(1, p))
  n_alphas <- length(alphas)

  # --- Subsampling (parallelised) ------------------------------------------
  raw <- array(
    unlist(
      mclapply(
        seq_len(B),
        function(b) .subsample(X, y, alphas),
        mc.cores = n_jobs
      )
    ),
    dim = c(n_alphas, 2L, p, B)
  )

  # Reshape to (n_alphas, 2B, p)
  Z <- array(NA_real_, dim = c(n_alphas, 2L * B, p))
  for (b in seq_len(B)) {
    Z[, (2L * (b - 1L) + 1L):(2L * b), ] <- raw[, , , b]
  }

  # --- Drop alphas with too many non-convergences --------------------------
  converge_prop <- colMeans(
    sapply(seq_len(n_alphas),
           function(i) sapply(seq_len(2L * B),
                              function(j) as.integer(anyNA(Z[i, j, ]))))
  )
  converge <- converge_prop < 0.45

  if (any(!converge)) {
    max_index <- min(which(!converge)) - 1L
    Z        <- Z[seq_len(max_index), , , drop = FALSE]
    n_alphas <- max_index
    alphas   <- alphas[seq_len(max_index)]
  }

  # --- Stability paths -----------------------------------------------------
  average_selected <- sapply(
    seq_len(n_alphas),
    function(i) mean(rowSums(Z[i, , ], na.rm = TRUE))
  )

  stability_paths <- matrix(NA_real_, nrow = n_alphas, ncol = p)
  for (i in seq_len(n_alphas)) {
    stability_paths[i, ] <- colMeans(Z[i, , ], na.rm = TRUE)
  }

  # --- IPSS scores ---------------------------------------------------------
  scores_info <- .ipss_scores(
    stability_paths, B, alphas, average_selected, ipss_function, delta, cutoff
  )

  efp_scores_raw <- round(
    scores_info$integral / pmax(scores_info$scores, scores_info$integral / p),
    8L
  )
  efp_scores <- as.list(setNames(efp_scores_raw, seq_len(p)))
  q_values   <- .compute_qvalues(efp_scores)

  # --- Feature selection ---------------------------------------------------
  if (!is.null(target_fdr)) {
    selected <- as.integer(names(Filter(function(q) q <= target_fdr, q_values)))
  } else if (!is.null(target_fp)) {
    selected <- as.integer(names(Filter(function(s) s <= target_fp, efp_scores)))
  } else {
    selected <- integer(0)
  }

  new_ipss_result(
    selected_features = selected,
    stability_paths   = stability_paths,
    efp_scores        = efp_scores,
    q_values          = q_values
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

#' Compute IPSS scores along the regularisation path
#'
#' @param stability_paths Matrix (n_alphas x p) of selection frequencies.
#' @param B Integer. Number of subsampling iterations.
#' @param alphas Numeric vector of regularisation values.
#' @param average_selected Numeric vector of average selected features per
#'   alpha.
#' @param ipss_function Character. One of \code{"h1"}, \code{"h2"},
#'   \code{"h3"}.
#' @param delta Numeric. Integration weight exponent.
#' @param cutoff Numeric. Early-stop threshold for the integrated bound.
#'
#' @return A list with \code{scores} (numeric vector, length p),
#'   \code{integral} (numeric scalar), and \code{stop_index} (integer).
#'
#' @keywords internal
.ipss_scores <- function(stability_paths, B, alphas, average_selected,
                         ipss_function, delta, cutoff) {

  if (!ipss_function %in% c("h1", "h2", "h3")) {
    stop(sprintf("ipss_function must be 'h1', 'h2', or 'h3', got '%s'.",
                 ipss_function))
  }

  p <- ncol(stability_paths)
  m <- switch(ipss_function, h1 = 1L, h2 = 2L, h3 = 3L)

  h_m <- function(x) ifelse(x <= 0.5, 0, (2 * x - 1)^m)

  integral_result <- switch(
    as.character(m),
    "1" = .integrate_f(average_selected^2 / p, alphas, delta, cutoff),
    "2" = {
      .integrate_f(average_selected^2 / (p * B) +
                     (B - 1) * average_selected^4 / (B * p^3),
                   alphas, delta, cutoff)
    },
    "3" = {
      .integrate_f(average_selected^2 / (p * B^2) +
                     (3 * (B - 1) * average_selected^4) / (p^3 * B^2) +
                     ((B - 1) * (B - 2) * average_selected^6) / (p^5 * B^2),
                   alphas, delta, cutoff)
    }
  )

  integral    <- integral_result$output
  stop_index  <- integral_result$stop_index
  alphas_stop <- alphas[seq_len(stop_index)]

  scores <- vapply(seq_len(p), function(i) {
    .integrate_f(h_m(stability_paths[seq_len(stop_index), i]),
                 alphas_stop, delta)$output
  }, numeric(1))

  list(scores = scores, integral = integral, stop_index = stop_index)
}


#' One subsampling iteration
#'
#' Splits data in half and fits LASSO Fine-Gray on each half.
#'
#' @param X Numeric matrix (n x p).
#' @param y Two-column matrix (time, status).
#' @param alphas Numeric vector of regularisation values.
#'
#' @return Array (n_alphas x 2 x p) of binary selection indicators.
#'   \code{NA} for non-converged fits.
#'
#' @keywords internal
.subsample <- function(X, y, alphas) {

  n          <- nrow(X)
  p          <- ncol(X)
  indices    <- sample.int(n)
  n_split    <- floor(n / 2L)
  indicators <- array(NA_real_, dim = c(length(alphas), 2L, p))

  for (half in 0:1) {
    idx <- if (half == 0L) indices[seq_len(n_split)] else indices[(n_split + 1L):n]
    indicators[, half + 1L, ] <- .fit_lasso_fg(X[idx, ], y[idx, ], alphas)
  }

  indicators
}


#' Fit LASSO Fine-Gray and return binary selection matrix
#'
#' @param X Numeric matrix (n_sub x p).
#' @param y Two-column matrix (time, status).
#' @param alphas Numeric vector of lambda values.
#'
#' @return Integer matrix (n_alphas x p): 1 = selected, 0 = not, NA =
#'   non-convergence.
#'
#' @keywords internal
.fit_lasso_fg <- function(X, y, alphas) {

  survobj <- fastcmprsk::Crisk(ftime = y[, 1], fstatus = y[, 2])

  model <- fastcmprsk::fastCrrp(
    survobj ~ .,
    data            = data.frame(X),
    penalty         = "LASSO",
    lambda          = alphas,
    alpha           = 0.95,
    penalty.factor  = rep(1, ncol(X)),
    standardize     = FALSE,
    getBreslowJumps = FALSE
  )

  coefs <- t(model$coef)
  coefs[!as.logical(model$converged), ] <- NA_real_
  matrix(as.integer(coefs != 0), nrow = nrow(coefs), ncol = ncol(coefs))
}


#' Compute the regularisation path
#'
#' Finds lambda_max via bisection and builds a log-spaced grid. Truncates at
#' the first lambda where the model fails to converge.
#'
#' @param X Numeric matrix (n x p).
#' @param y Two-column matrix (time, status).
#' @param n_alphas Integer. Desired path length.
#' @param penalty.factor Numeric vector (length p).
#'
#' @return Numeric vector of lambda values.
#'
#' @keywords internal
.compute_alphas <- function(X, y, n_alphas, penalty.factor) {

  ratio    <- 0.05
  max_iter <- 3L
  iter     <- 0L
  path_ok  <- FALSE
  alphas   <- NULL

  alpha_max <- .get_alpha_max(X, y, penalty.factor)
  survobj   <- fastcmprsk::Crisk(ftime = y[, 1], fstatus = y[, 2])

  while (!path_ok && iter < max_iter) {

    alpha_min <- ratio * alpha_max
    alphas    <- alpha_max *
      (alpha_min / alpha_max)^(seq(0, n_alphas - 1L) / (n_alphas - 1L))

    for (i in rev(seq_len(n_alphas - 3L))) {
      model <- fastcmprsk::fastCrrp(
        survobj ~ .,
        data            = data.frame(X),
        penalty         = "LASSO",
        lambda          = alphas[i],
        penalty.factor  = penalty.factor,
        standardize     = FALSE,
        getBreslowJumps = FALSE
      )
      if (!as.logical(model$converged)) {
        alphas <- alphas[seq_len(i)]
        break
      }
    }

    if (length(alphas) < n_alphas) {
      ratio <- ratio * 2
      iter  <- iter + 1L
    } else {
      path_ok <- TRUE
    }
  }

  alphas
}


#' Binary search for lambda_max
#'
#' Finds the smallest lambda such that no feature is selected.
#'
#' @param X Numeric matrix (n x p).
#' @param y Two-column matrix (time, status).
#' @param penalty.factor Numeric vector (length p).
#' @param eps Numeric. Bisection tolerance. Default \code{0.01}.
#' @param iter_max Integer. Maximum iterations. Default \code{40}.
#'
#' @return Numeric scalar: lambda_max.
#'
#' @keywords internal
.get_alpha_max <- function(X, y, penalty.factor, eps = 0.01, iter_max = 40L) {

  alpha_sup <- 50
  alpha_inf <- 0
  survobj   <- fastcmprsk::Crisk(ftime = y[, 1], fstatus = y[, 2])
  iter      <- 0L

  repeat {
    alpha <- (alpha_sup + alpha_inf) / 2
    model <- fastcmprsk::fastCrrp(
      survobj ~ .,
      data            = data.frame(X),
      penalty         = "LASSO",
      lambda          = alpha,
      alpha           = 0.95,
      standardize     = FALSE,
      penalty.factor  = penalty.factor,
      getBreslowJumps = FALSE
    )
    if (sum(model$coef != 0) == 0L) alpha_sup <- alpha else alpha_inf <- alpha
    iter <- iter + 1L
    if ((alpha_sup - alpha_inf) <= eps || iter >= iter_max) break
  }

  alpha_sup
}


#' Compute q-values from EFP scores
#'
#' For each feature, the q-value is the minimum FDR achievable when including
#' that feature in the selected set.
#'
#' @param efp_scores Named list of EFP scores (one per feature).
#'
#' @return Named numeric vector of q-values.
#'
#' @keywords internal
.compute_qvalues <- function(efp_scores) {

  efp_vec <- unlist(efp_scores)

  fdrs <- lapply(efp_vec, function(t) {
    subset_leq <- efp_vec[efp_vec <= t]
    c(threshold = t, fdr = max(subset_leq) / length(subset_leq))
  })

  vapply(names(efp_scores), function(feat) {
    score <- efp_scores[[feat]]
    min(vapply(fdrs,
               function(x) if (score <= x["threshold"]) x["fdr"] else Inf,
               numeric(1)))
  }, numeric(1))
}


#' Numerical integration along the regularisation path
#'
#' Weighted integral of \code{values} w.r.t. \code{alphas} on a log scale.
#' Stops early if \code{cutoff} is exceeded.
#'
#' @param values Numeric vector (length n_alphas).
#' @param alphas Numeric vector of lambda values (decreasing).
#' @param delta Numeric. Weight exponent (\code{1} = uniform on log scale).
#' @param cutoff Numeric or \code{NULL}. Early-stop threshold.
#'
#' @return List with \code{output} (numeric scalar) and \code{stop_index}
#'   (integer).
#'
#' @keywords internal
.integrate_f <- function(values, alphas, delta = 1, cutoff = NULL) {

  n_alphas <- length(alphas)
  a <- min(alphas)
  b <- max(alphas)

  normalization <- if (delta == 1) {
    (1 - (a / b)^(1 / n_alphas)) / log(b / a)
  } else {
    (1 - delta) * (1 - (a / b)^(1 / n_alphas)) / (b^(1 - delta) - a^(1 - delta))
  }

  output     <- 0
  stop_index <- n_alphas

  for (i in seq(2L, n_alphas)) {
    weight         <- if (delta == 1) 1 else alphas[i]^(1 - delta)
    updated_output <- output + normalization * weight * values[i - 1L]
    if (!is.null(cutoff) && updated_output > cutoff) {
      stop_index <- i
      break
    }
    output <- updated_output
  }

  list(output = output, stop_index = stop_index - 1L)
}
