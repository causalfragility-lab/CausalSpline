#' CausalSpline: Nonlinear Causal Dose-Response Estimation via Splines
#'
#' @description
#' CausalSpline estimates causal dose-response functions
#' \eqn{E[Y(t)]} for continuous treatments under unconfoundedness, without
#' imposing linearity on the treatment effect. The package fills a gap in the
#' causal inference ecosystem: most tools assume a scalar \eqn{\beta_1 T}
#' treatment effect, whereas real policy or health applications frequently
#' exhibit thresholds, diminishing returns, and non-monotone relationships.
#'
#' @section Core functions:
#' \describe{
#'   \item{\code{\link{causal_spline}}}{Main fitting function. Supports IPW,
#'     G-computation, and doubly-robust estimation.}
#'   \item{\code{\link{gps_model}}}{Fit the generalised propensity score model
#'     for continuous treatments.}
#'   \item{\code{\link{trim_weights}}}{Winsorise extreme IPW weights.}
#'   \item{\code{\link{check_overlap}}}{Diagnose positivity (ESS, weight plots).}
#'   \item{\code{\link{dose_response_curve}}}{Extract the curve data frame.}
#'   \item{\code{\link{simulate_dose_response}}}{Simulate nonlinear dose-response
#'     data for benchmarking.}
#'   \item{\code{\link{gradient_curve}}}{Numerical first and second derivatives
#'     of the fitted dose-response curve.}
#'   \item{\code{\link{fragility_curve}}}{Pointwise fragility with adaptive eps,
#'     interpretation zones, and uncertainty normalisation.}
#'   \item{\code{\link{region_fragility}}}{Integrated fragility over a treatment
#'     interval.}
#' }
#'
#' @section Causal identification:
#' Identification relies on two standard assumptions:
#' \enumerate{
#'   \item \strong{Unconfoundedness} (strong ignorability):
#'     \eqn{Y(t) \perp T \mid X} for all \eqn{t}.
#'   \item \strong{Positivity} (overlap):
#'     \eqn{f(t \mid X = x) > 0} for all \eqn{t} in the support of \eqn{T}
#'     and all \eqn{x} in the support of \eqn{X}.
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{IPW}{The outcome is regressed on a spline of T using GPS-based
#'     inverse probability weights. Consistent if the GPS model is correctly
#'     specified.}
#'   \item{G-computation}{The outcome model includes spline(T) + X. The
#'     curve is obtained by marginalising over the observed X distribution.
#'     Consistent if the outcome model is correctly specified.}
#'   \item{Doubly Robust}{Combines IPW and g-computation. Consistent if
#'     at least one of the two models is correctly specified.}
#' }
#'
#' @references
#' Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous
#' treatments. \doi{10.1002/0470090456.ch7}
#'
#' Imbens, G. W. (2000). The role of the propensity score in estimating
#' dose-response functions. \emph{Biometrika}, 87(3), 706-710.
#' \doi{10.1093/biomet/87.3.706}
#'
#' @importFrom sandwich vcovHC vcovCL
#' @importFrom boot boot boot.ci
#' @importFrom splines ns bs
#' @importFrom stats lm glm coef predict fitted residuals sd var quantile
#'   qnorm dnorm dgamma dbeta rnorm rbinom Gamma quasibinomial
#'   model.frame model.matrix model.response as.formula complete.cases
#'   na.omit weighted.mean setNames median approx density
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_rug geom_hline
#'   geom_histogram geom_area annotate labs theme_bw theme element_text
#'   element_blank
#' @importFrom utils head tail
"_PACKAGE"

utils::globalVariables(c("t", "estimate", "lower", "upper",
                         "true_effect", "w", "T", "fragility",
                         "high_fragility", "dens_scaled",
                         "support_density", "fragility_type"))
