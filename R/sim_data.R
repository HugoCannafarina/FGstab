#' Simulate competing risks data under the Fine-Gray model
#'
#' Generates a synthetic dataset for a competing risks model with two event
#' types (event of interest and competing event), following the
#' sub-distribution hazard model of Fine and Gray (1999). Covariates
#' \eqn{X \in \mathbb{R}^{n \times p}} are drawn from a multivariate normal
#' distribution \eqn{X \sim \mathcal{N}(0, \Sigma)} whose covariance structure
#' \eqn{\Sigma} is controlled by the \code{design} argument.
#'
#' @section Covariance structures:
#'
#' \describe{
#'   \item{\strong{Independent} (\code{"independent"})}{
#'     \eqn{\Sigma = I_p}.
#'     Each covariate is drawn i.i.d. from \eqn{\mathcal{N}(0, 1)}.
#'   }
#'   \item{\strong{Toeplitz} (\code{"toeplitz"})}{
#'     \deqn{\Sigma_{jk} = \rho^{|j - k|}, \quad j, k = 1, \ldots, p.}
#'     Correlations decay exponentially with the distance between indices.
#'     The parameter \code{rho} \eqn{\in (0, 1)} controls the decay rate:
#'     values close to 1 yield strongly correlated neighbouring covariates,
#'     values close to 0 approach the independent design.
#'   }
#'   \item{\strong{Block} (\code{"block"})}{
#'     \deqn{
#'       \Sigma_{jk} = \begin{cases}
#'         1   & \text{if } j = k, \\
#'         0.5 & \text{if } j \neq k \text{ and }
#'                        \lfloor j / B \rfloor = \lfloor k / B \rfloor, \\
#'         0   & \text{otherwise,}
#'       \end{cases}
#'     }
#'     where \eqn{B = 10} is the block size. Covariates within the same block
#'     share a fixed pairwise correlation of 0.5; covariates in different
#'     blocks are independent.
#'   }
#' }
#'
#' @section Data-generating mechanism:
#'
#' True coefficients \eqn{\beta \in \mathbb{R}^p} are sparse: \code{n_pred}
#' indices are drawn uniformly at random and assigned coefficients from
#' \eqn{\mathcal{U}(-1, 1)}; all other coefficients are zero. The marginal
#' probability of the event of interest for individual \eqn{i} is
#' \deqn{
#'   \pi_i = 1 - (1 - \pi_0)^{\exp(X_i^\top \beta)},
#' }
#' where \eqn{\pi_0 = 0.75} is the baseline event rate. Event times for the
#' event of interest are drawn via quantile inversion of the sub-distribution
#' function; competing event times are drawn from
#' \eqn{\mathrm{Exp}(\text{rate} = 3)}. Censoring is applied by thresholding
#' at the \code{censor_rate} quantile of the observed times.
#'
#' @param n_ind Integer. Number of individuals. Default \code{250}.
#' @param n_pred Integer. Number of truly predictive covariates
#'   (\eqn{|\mathcal{S}|}). Default \code{10}.
#' @param n_var Integer. Total number of covariates (\eqn{p}). Default
#'   \code{5000}.
#' @param censor_rate Numeric in \eqn{[0, 1)}. Target proportion of censored
#'   observations. Default \code{0.3}.
#' @param design Character string. Covariance structure of \eqn{\Sigma}.
#'   One of \code{"independent"}, \code{"toeplitz"}, or \code{"block"}.
#'   See the \emph{Covariance structures} section for details.
#' @param rho Numeric in \eqn{(0, 1)}. Decay parameter for the
#'   \code{"toeplitz"} design. Ignored for other designs. Default \code{0.9}.
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{\code{X}}{Numeric matrix of covariates (\code{n_ind} x
#'       \code{n_var}), centred and scaled.}
#'     \item{\code{surv}}{Data frame with columns \code{Tobs} (observed time),
#'       \code{evtype} (0 = censored, 1 = event of interest, 2 = competing
#'       event), and \code{id}.}
#'     \item{\code{true_pred_index}}{Integer vector of length \code{n_pred}:
#'       indices of the truly predictive covariates (\eqn{\mathcal{S}}).}
#'     \item{\code{beta}}{Numeric vector of true coefficients (length
#'       \code{n_var}; non-zero only at \code{true_pred_index}).}
#'   }
#'
#' @references
#' Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the
#' subdistribution of a competing risk. \emph{Journal of the American
#' Statistical Association}, \strong{94}(446), 496--509.
#'
#' @importFrom MASS mvrnorm
#'
#' @examples
#' # Independent covariates
#' dat <- sim_competing(n_ind = 100, n_pred = 5, n_var = 200, design = "independent")
#' str(dat)
#'
#' # Toeplitz covariance (moderate correlation)
#' dat_toep <- sim_competing(n_ind = 150, n_pred = 8, n_var = 300,
#'                      design = "toeplitz", rho = 0.7)
#'
#' # Block covariance
#' dat_block <- sim_competing(n_ind = 150, n_pred = 8, n_var = 300,
#'                       design = "block")
#'
#' @export
sim_competing <- function(n_ind       = 250L,
                     n_pred      = 10L,
                     n_var       = 5000L,
                     censor_rate = 0.3,
                     design      = c("independent", "toeplitz", "block"),
                     rho         = 0.9) {

  design        <- match.arg(design)
  main_evt_rate <- 0.75

  # --- Covariate matrix ----------------------------------------------------
  if (design == "independent") {

    X <- matrix(rnorm(n_var * n_ind), nrow = n_ind, ncol = n_var)

  } else if (design == "toeplitz") {

    idx   <- seq_len(n_var)
    Sigma <- rho^abs(outer(idx, idx, "-"))
    diag(Sigma) <- 1
    X <- MASS::mvrnorm(n = n_ind, mu = rep(0, n_var), Sigma = Sigma)

  } else {  # "block"

    n_block   <- 10L
    cov_value <- 0.5
    Sigma     <- matrix(0, nrow = n_var, ncol = n_var)
    for (k in seq_len(n_var)) {
      for (m in seq_len(n_var)) {
        if (k %% n_block == m %% n_block) Sigma[k, m] <- cov_value
      }
    }
    diag(Sigma) <- 1
    X <- MASS::mvrnorm(n_ind, mu = rep(0, n_var), Sigma = Sigma)

  }

  X <- scale(X)

  # --- True coefficients ---------------------------------------------------
  true_pred_index          <- sample(seq_len(n_var), n_pred)
  beta                     <- numeric(n_var)
  beta[true_pred_index]    <- runif(n_pred, -1, 1)

  # --- Event times ---------------------------------------------------------
  p_main   <- 1 - (1 - main_evt_rate)^exp(X %*% beta)
  evt_type <- rbinom(n_ind, size = 1L, prob = 1 - p_main) + 1L

  evt_time <- numeric(n_ind)
  evt_time[evt_type == 1L] <- .quantile_main(
    X[evt_type == 1L, , drop = FALSE],
    beta, p_main[evt_type == 1L, , drop = FALSE],
    main_evt_rate
  )
  evt_time[evt_type == 2L] <- rexp(sum(evt_type == 2L), rate = 3)

  # --- Censoring -----------------------------------------------------------
  censor_threshold         <- sort(evt_time, decreasing = TRUE)[
    max(1L, round(n_ind * censor_rate))
  ]
  evt_type[evt_time >= censor_threshold] <- 0L

  list(
    X               = X,
    surv            = data.frame(Tobs   = evt_time,
                                 evtype = evt_type,
                                 id     = seq_len(n_ind)),
    true_pred_index = true_pred_index,
    beta            = beta
  )
}


#' Simulate survival data under the Cox proportional hazards model
#'
#' Generates a synthetic dataset for a standard survival model with no
#' competing risks, following the Cox proportional hazards model with an
#' exponential baseline hazard. Covariates
#' \eqn{X \in \mathbb{R}^{n \times p}} are drawn from a multivariate normal
#' distribution \eqn{X \sim \mathcal{N}(0, \Sigma)} whose covariance structure
#' \eqn{\Sigma} is controlled by the \code{design} argument (same options as
#' \code{\link{sim_competing}}).
#'
#' @section Data-generating mechanism:
#'
#' True coefficients \eqn{\beta \in \mathbb{R}^p} are sparse: \code{n_pred}
#' indices are drawn uniformly at random and assigned coefficients from
#' \eqn{\mathcal{U}(-1, 1)}; all other coefficients are zero. Event times are
#' drawn from an exponential distribution with rate \eqn{\exp(X_i^\top \beta)}:
#' \deqn{T_i \sim \mathrm{Exp}\!\left(\exp(X_i^\top \beta)\right).}
#' This corresponds to a Cox model with unit exponential baseline hazard.
#' Censoring is applied by thresholding at the \code{censor_rate} quantile of
#' the observed times.
#'
#' The output is formatted identically to \code{\link{sim_competing}} with
#' \code{evtype} \eqn{\in \{0, 1\}} (no competing event), so it can be passed
#' directly to \code{\link{FGstab}}, which then operates as IPSS for the Cox
#' model.
#'
#' @inheritParams sim_competing
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{\code{X}}{Numeric matrix of covariates (\code{n_ind} x
#'       \code{n_var}), centred and scaled.}
#'     \item{\code{surv}}{Data frame with columns \code{Tobs} (observed time),
#'       \code{evtype} (0 = censored, 1 = event), and \code{id}.}
#'     \item{\code{true_pred_index}}{Integer vector of length \code{n_pred}:
#'       indices of the truly predictive covariates (\eqn{\mathcal{S}}).}
#'     \item{\code{beta}}{Numeric vector of true coefficients (length
#'       \code{n_var}; non-zero only at \code{true_pred_index}).}
#'   }
#'
#' @seealso \code{\link{sim_competing}}, \code{\link{FGstab}}
#'
#' @importFrom MASS mvrnorm
#'
#' @examples
#' dat <- sim_survival(n_ind = 100, n_pred = 5, n_var = 200,
#'                     design = "independent")
#' str(dat)
#'
#' \dontrun{
#' y <- as.matrix(dat$surv[, c("Tobs", "evtype")])
#' res <- FGstab(dat$X, y, n_alphas = 10, target_fdr = 0.2)
#' print(res)
#' }
#'
#' @export
sim_survival <- function(n_ind       = 250L,
                         n_pred      = 10L,
                         n_var       = 5000L,
                         censor_rate = 0.3,
                         design      = c("independent", "toeplitz", "block"),
                         rho         = 0.9) {

  design <- match.arg(design)

  # --- Covariate matrix ----------------------------------------------------
  if (design == "independent") {

    X <- matrix(rnorm(n_var * n_ind), nrow = n_ind, ncol = n_var)

  } else if (design == "toeplitz") {

    idx   <- seq_len(n_var)
    Sigma <- rho^abs(outer(idx, idx, "-"))
    diag(Sigma) <- 1
    X <- MASS::mvrnorm(n = n_ind, mu = rep(0, n_var), Sigma = Sigma)

  } else {  # "block"

    n_block   <- 10L
    cov_value <- 0.5
    Sigma     <- matrix(0, nrow = n_var, ncol = n_var)
    for (k in seq_len(n_var)) {
      for (m in seq_len(n_var)) {
        if (k %% n_block == m %% n_block) Sigma[k, m] <- cov_value
      }
    }
    diag(Sigma) <- 1
    X <- MASS::mvrnorm(n_ind, mu = rep(0, n_var), Sigma = Sigma)

  }

  X <- scale(X)

  # --- True coefficients ---------------------------------------------------
  true_pred_index          <- sample(seq_len(n_var), n_pred)
  beta                     <- numeric(n_var)
  beta[true_pred_index]    <- runif(n_pred, -1, 1)

  # --- Event times: Cox model with exponential baseline --------------------
  evt_time <- rexp(n_ind, rate = as.vector(exp(X %*% beta)))
  evt_type <- rep(1L, n_ind)

  # --- Censoring -----------------------------------------------------------
  censor_threshold         <- sort(evt_time, decreasing = TRUE)[
    max(1L, round(n_ind * censor_rate))
  ]
  evt_type[evt_time >= censor_threshold] <- 0L

  list(
    X               = X,
    surv            = data.frame(Tobs   = evt_time,
                                 evtype = evt_type,
                                 id     = seq_len(n_ind)),
    true_pred_index = true_pred_index,
    beta            = beta
  )
}


# =============================================================================
# Internal helper
# =============================================================================

#' Simulate event times for the event of interest via quantile inversion
#'
#' Draws event times from the sub-distribution of the event of interest using
#' the inverse-CDF method under the Fine-Gray model.
#'
#' @param X_main Numeric matrix of covariates for individuals with the event
#'   of interest.
#' @param beta Numeric vector of true coefficients (length p).
#' @param p_main Numeric vector of marginal probabilities of the main event
#'   for the same individuals.
#' @param main_evt_rate Numeric. Overall rate of the main event in the
#'   population.
#'
#' @return Numeric vector of simulated event times.
#'
#' @keywords internal
.quantile_main <- function(X_main, beta, p_main, main_evt_rate) {
  u <- runif(nrow(X_main))
  -log(1 - (1 - (1 - u * p_main)^exp(-X_main %*% beta)) / main_evt_rate)
}
