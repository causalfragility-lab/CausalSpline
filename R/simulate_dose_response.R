#' Simulate nonlinear dose-response data
#'
#' @description
#' Generates synthetic observational data with a continuous treatment,
#' confounders, and a nonlinear dose-response outcome. Useful for testing,
#' benchmarking, and illustrating the CausalSpline package.
#'
#' @param n Integer. Sample size. Default \code{500}.
#' @param dgp Character string. Data-generating process:
#' \describe{
#'   \item{\code{"threshold"}}{Flat below treatment = 3, steep rise above.}
#'   \item{\code{"diminishing"}}{Concave relationship with diminishing returns.}
#'   \item{\code{"nonmonotone"}}{Inverted-U relationship.}
#'   \item{\code{"linear"}}{Standard linear effect (useful as baseline).}
#'   \item{\code{"sinusoidal"}}{Oscillatory effect (difficult, high df needed).}
#' }
#' @param confounding Numeric scalar in \code{[0, 1]}. Strength of confounding
#'   (correlation between treatment and confounders). Default \code{0.5}.
#' @param sigma_y Numeric. Standard deviation of the outcome noise.
#'   Default \code{1}.
#' @param seed Integer or \code{NULL}. Random seed. Default \code{NULL}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{Y}}{Observed outcome.}
#'   \item{\code{T}}{Observed (confounded) treatment.}
#'   \item{\code{X1}, \code{X2}, \code{X3}}{Confounders.}
#'   \item{\code{Y0}}{Potential outcome at T = 0 (for evaluation).}
#'   \item{\code{true_effect}}{f(T) at each observed T value.}
#' }
#'
#' @examples
#' dat <- simulate_dose_response(n = 500, dgp = "threshold", seed = 1)
#' plot(dat$T, dat$true_effect, type = "l",
#'      xlab = "Treatment", ylab = "True causal effect")
#'
#' dat2 <- simulate_dose_response(n = 1000, dgp = "nonmonotone",
#'                                 confounding = 0.8, seed = 42)
#' hist(dat2$T, main = "Treatment distribution", xlab = "T")
#'
#' @export
simulate_dose_response <- function(n           = 500L,
                                   dgp         = c("threshold", "diminishing",
                                                   "nonmonotone", "linear",
                                                   "sinusoidal"),
                                   confounding = 0.5,
                                   sigma_y     = 1,
                                   seed        = NULL) {
  dgp <- match.arg(dgp)
  if (!is.null(seed)) set.seed(seed)

  # ── Confounders ──────────────────────────────────────────────────────────────
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rbinom(n, 1, 0.5)

  # ── Treatment (continuous, confounded by X) ──────────────────────────────────
  # Treatment range approximately [0, 10]
  T_star <- 5 + confounding * (0.8 * X1 - 0.5 * X2 + 0.4 * X3) +
            sqrt(1 - confounding^2) * rnorm(n)
  Tr <- pmax(pmin(T_star, 10), 0)   # soft clip to [0, 10]

  # ── True dose-response function f(T) ─────────────────────────────────────────
  f_T <- switch(dgp,
    threshold   = ifelse(Tr < 3, 0, 2 * (Tr - 3)),
    diminishing = 5 * log(1 + Tr),
    nonmonotone = -0.2 * (Tr - 5)^2 + 5,
    linear      = 0.8 * Tr,
    sinusoidal  = 3 * sin(Tr * pi / 5) + 0.3 * Tr
  )

  # ── Outcome (confounded via direct X effect) ─────────────────────────────────
  Y <- f_T + 0.5 * X1 - 0.3 * X2 + 0.2 * X3 + rnorm(n, sd = sigma_y)

  data.frame(
    Y           = Y,
    T           = Tr,
    X1          = X1,
    X2          = X2,
    X3          = X3,
    true_effect = f_T
  )
}
