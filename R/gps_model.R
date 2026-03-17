#' Fit the Generalised Propensity Score model
#'
#' @description
#' Fits a model for the conditional distribution of a continuous treatment
#' \eqn{T | X} (the generalised propensity score). Covariates are entered
#' linearly by default; spline transformations of the treatment are not needed
#' here since this models the treatment, not the outcome.
#'
#' @param treatment Numeric vector of treatment values.
#' @param covariates Numeric matrix of confounders.
#' @param spline_type Character. \code{"ns"} or \code{"bs"} for spline
#'   transformations of continuous covariates. Default \code{"ns"}.
#' @param df Integer. Degrees of freedom for covariate splines. Default
#'   \code{4}. Set to \code{NULL} to use linear covariates only.
#' @param family Character. Distribution for the GPS model.
#'   \code{"gaussian"} (default), \code{"gamma"}, or \code{"beta"}.
#' @param verbose Logical. Default \code{FALSE}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{model}}{Fitted model object (\code{lm} or \code{glm}).}
#'   \item{\code{family}}{Family used.}
#'   \item{\code{r_squared}}{R-squared of the GPS model (diagnostic).}
#' }
#'
#' @examples
#' dat <- simulate_dose_response(n = 400, dgp = "linear")
#' X   <- as.matrix(dat[, c("X1", "X2", "X3")])
#' gps <- gps_model(dat$T, X)
#' summary(gps$model)
#'
#' @export
gps_model <- function(treatment,
                      covariates,
                      spline_type = c("ns", "bs"),
                      df          = 4L,
                      family      = c("gaussian", "gamma", "beta"),
                      verbose     = FALSE) {

  spline_type <- match.arg(spline_type)
  family      <- match.arg(family)

  X_df  <- as.data.frame(covariates)
  df_gps <- data.frame(T_treat = treatment, X_df)

  if (family == "gaussian") {
    fit <- lm(T_treat ~ ., data = df_gps)
    r2  <- summary(fit)$r.squared
  } else if (family == "gamma") {
    fit <- glm(T_treat ~ ., data = df_gps, family = Gamma(link = "log"))
    r2  <- 1 - fit$deviance / fit$null.deviance
  } else {
    # beta: rescale treatment to (0, 1)
    tr_sc <- (treatment - min(treatment) + 1e-6) /
             (max(treatment) - min(treatment) + 2e-6)
    df_gps$T_treat <- tr_sc
    fit <- glm(T_treat ~ ., data = df_gps,
               family = quasibinomial(link = "logit"))
    r2  <- 1 - fit$deviance / fit$null.deviance
  }

  if (verbose) {
    message("[GPS] family = ", family,
            " | R^2 = ", round(r2, 3))
  }

  list(model = fit, family = family, r_squared = r2)
}
