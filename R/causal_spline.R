#' Nonlinear causal dose-response estimation via splines
#'
#' @description
#' Estimates the causal dose-response function \eqn{E[Y(t)]} for a continuous
#' treatment \eqn{T} under unconfoundedness, using either Inverse Probability
#' Weighting (IPW) with the generalised propensity score (GPS) or
#' G-computation (outcome regression). The exposure-response curve is modelled
#' as a natural cubic spline or B-spline, allowing thresholds, diminishing
#' returns, and other nonlinear patterns to be recovered without parametric
#' assumptions on the functional form.
#'
#' @param formula A two-sided formula of the form \code{outcome ~ treatment | covariate1 + covariate2 + ...}.
#'   The vertical bar \code{|} separates the treatment variable from the
#'   confounders. Example: \code{Y ~ T | X1 + X2 + X3}.
#' @param data A data frame containing all variables in \code{formula}.
#' @param method Character string. Estimation method: \code{"ipw"} (inverse
#'   probability weighting via GPS), \code{"gcomp"} (g-computation / outcome
#'   regression), or \code{"dr"} (doubly-robust, combines both). Default
#'   \code{"ipw"}.
#' @param spline_type Character string. Type of spline basis for \eqn{f(T)}:
#'   \code{"ns"} (natural cubic spline, default) or \code{"bs"} (B-spline).
#' @param df_exposure Integer. Degrees of freedom for the treatment spline
#'   \eqn{f(T)}. Default \code{4}.
#' @param knots_exposure Numeric vector of interior knot positions for the
#'   treatment spline. If \code{NULL} (default), knots are placed at quantiles
#'   of the treatment distribution.
#' @param df_gps Integer. Degrees of freedom for the GPS (propensity) spline
#'   used in the \code{"ipw"} method. Default \code{4}.
#' @param gps_family Character string. Family for the GPS model:
#'   \code{"gaussian"} (default) for a linear GPS, \code{"gamma"}, or
#'   \code{"beta"} for positive / bounded treatments.
#' @param trim_quantiles Numeric vector of length 2 giving lower and upper
#'   quantiles for weight trimming. Default \code{c(0.01, 0.99)}. Set to
#'   \code{NULL} to skip trimming.
#' @param eval_grid Integer. Number of equally-spaced treatment values at
#'   which to evaluate \eqn{E[Y(t)]}. Default \code{100}.
#' @param eval_points Numeric vector of specific treatment values at which to
#'   evaluate the curve. Overrides \code{eval_grid} if supplied.
#' @param se_method Character string. Method for standard errors:
#'   \code{"sandwich"} (default, fast) or \code{"bootstrap"}.
#' @param boot_reps Integer. Number of bootstrap replications when
#'   \code{se_method = "bootstrap"}. Default \code{500}.
#' @param conf_level Numeric. Confidence level for intervals. Default
#'   \code{0.95}.
#' @param verbose Logical. Print progress messages. Default \code{FALSE}.
#'
#' @return An object of class \code{"causal_spline"} with elements:
#' \describe{
#'   \item{\code{curve}}{A data frame with columns \code{t} (treatment grid),
#'     \code{estimate} (E[Y(t)]), \code{se}, \code{lower}, \code{upper}.}
#'   \item{\code{ate}}{Estimated average treatment effect over the observed
#'     treatment range (scalar).}
#'   \item{\code{weights}}{Numeric vector of final (trimmed, normalised)
#'     IPW weights (only for \code{method = "ipw"} or \code{"dr"}).}
#'   \item{\code{gps_fit}}{The fitted GPS model object.}
#'   \item{\code{outcome_fit}}{The fitted outcome model object.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{method}}{The estimation method used.}
#'   \item{\code{spline_type}}{Spline type used.}
#'   \item{\code{df_exposure}}{Degrees of freedom for the exposure spline.}
#'   \item{\code{data_summary}}{Summary statistics of the treatment variable.}
#' }
#'
#' @references
#' Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous
#' treatments. \emph{Applied Bayesian Modeling and Causal Inference from
#' Incomplete-Data Perspectives}, 226164, 73-84.
#' \doi{10.1002/0470090456.ch7}
#'
#' Flores, C. A., Flores-Lagunes, A., Gonzalez, A., & Neumann, T. C. (2012).
#' Estimating the effects of length of exposure to instruction in a training
#' program: the case of job corps. \emph{Review of Economics and Statistics},
#' 94(1), 153-171. \doi{10.1162/REST_a_00177}
#'
#' Imbens, G. W. (2000). The role of the propensity score in estimating
#' dose-response functions. \emph{Biometrika}, 87(3), 706-710.
#' \doi{10.1093/biomet/87.3.706}
#'
#' @examples
#' # Simulate nonlinear dose-response data
#' set.seed(42)
#' dat <- simulate_dose_response(n = 500, dgp = "threshold")
#'
#' # Fit with IPW
#' fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat, method = "ipw",
#'                      df_exposure = 5)
#' summary(fit)
#' plot(fit)
#'
#' # Fit with G-computation
#' fit_gc <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat,
#'                         method = "gcomp", df_exposure = 5)
#' plot(fit_gc, truth = dat)
#'
#' @export
causal_spline <- function(formula,
                          data,
                          method         = c("ipw", "gcomp", "dr"),
                          spline_type    = c("ns", "bs"),
                          df_exposure    = 4L,
                          knots_exposure = NULL,
                          df_gps         = 4L,
                          gps_family     = c("gaussian", "gamma", "beta"),
                          trim_quantiles = c(0.01, 0.99),
                          eval_grid      = 100L,
                          eval_points    = NULL,
                          se_method      = c("sandwich", "bootstrap"),
                          boot_reps      = 500L,
                          conf_level     = 0.95,
                          verbose        = FALSE) {

  cl          <- match.call()
  method      <- match.arg(method)
  spline_type <- match.arg(spline_type)
  gps_family  <- match.arg(gps_family)
  se_method   <- match.arg(se_method)

  # ── 1. Parse formula ────────────────────────────────────────────────────────
  parsed <- .parse_causal_formula(formula, data)
  Y  <- parsed$Y
  Tr <- parsed$T
  X  <- parsed$X
  n  <- length(Y)

  if (verbose) message("[CausalSpline] n = ", n, " | method = ", method)

  # ── 2. Build treatment evaluation grid ──────────────────────────────────────
  t_range <- range(Tr, na.rm = TRUE)
  t_grid  <- if (!is.null(eval_points)) {
    sort(eval_points)
  } else {
    seq(t_range[1], t_range[2], length.out = eval_grid)
  }

  # ── 3. Build spline basis for treatment ─────────────────────────────────────
  spline_fn <- if (spline_type == "ns") splines::ns else splines::bs
  if (is.null(knots_exposure)) {
    knot_probs     <- seq(0, 1, length.out = df_exposure + 1)[-c(1, df_exposure + 1)]
    knots_exposure <- stats::quantile(Tr, probs = knot_probs, na.rm = TRUE)
  }
  T_basis      <- spline_fn(Tr, knots = knots_exposure,
                            Boundary.knots = t_range, intercept = FALSE)
  T_basis_grid <- stats::predict(T_basis, newx = t_grid)

  # ── 4. Fit GPS model (needed for IPW / DR) ──────────────────────────────────
  gps_fit <- NULL
  weights  <- rep(1, n)

  if (method %in% c("ipw", "dr")) {
    gps_fit <- gps_model(
      treatment   = Tr,
      covariates  = X,
      spline_type = spline_type,
      df          = df_gps,
      family      = gps_family,
      verbose     = verbose
    )
    weights <- .compute_ipw_weights(Tr, gps_fit, gps_family)
    if (!is.null(trim_quantiles)) {
      weights <- trim_weights(weights, trim_quantiles)
    }
    weights <- weights / mean(weights)  # Hajek-style normalisation
    if (verbose) message("[CausalSpline] ESS = ",
                         round(sum(weights)^2 / sum(weights^2)))
  }

  # ── 5. Estimate outcome model / g-computation ────────────────────────────────
  outcome_fit <- NULL
  curve_est   <- NULL

  if (method == "ipw") {
    curve_est   <- .fit_ipw_spline(Y, T_basis, T_basis_grid, weights,
                                   se_method, boot_reps, conf_level, Tr, t_grid)
    outcome_fit <- curve_est$fit

  } else if (method == "gcomp") {
    curve_est   <- .fit_gcomp_spline(Y, Tr, T_basis, T_basis_grid, X,
                                     spline_type, knots_exposure, t_range,
                                     t_grid, se_method, boot_reps, conf_level)
    outcome_fit <- curve_est$fit

  } else {  # doubly-robust
    ipw_part    <- .fit_ipw_spline(Y, T_basis, T_basis_grid, weights,
                                   se_method, boot_reps, conf_level, Tr, t_grid)
    gcomp_part  <- .fit_gcomp_spline(Y, Tr, T_basis, T_basis_grid, X,
                                     spline_type, knots_exposure, t_range,
                                     t_grid, se_method, boot_reps, conf_level)
    curve_est   <- .combine_dr(ipw_part, gcomp_part, conf_level)
    outcome_fit <- gcomp_part$fit
  }

  # ── 6. Assemble output ───────────────────────────────────────────────────────
  structure(
    list(
      curve        = curve_est$curve,
      ate          = mean(curve_est$curve$estimate),
      weights      = if (method != "gcomp") weights else NULL,
      gps_fit      = gps_fit,
      outcome_fit  = outcome_fit,
      call         = cl,
      method       = method,
      spline_type  = spline_type,
      df_exposure  = df_exposure,
      knots        = knots_exposure,
      t_range      = t_range,
      conf_level   = conf_level,
      data_summary = list(
        n      = n,
        t_mean = mean(Tr),
        t_sd   = stats::sd(Tr),
        t_range = t_range,
        y_mean = mean(Y)
      )
    ),
    class = "causal_spline"
  )
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.parse_causal_formula <- function(formula, data) {
  f_char <- deparse(formula)
  if (!grepl("\\|", f_char))
    stop("Formula must use '|' to separate treatment from covariates.\n",
         "  Example: Y ~ T | X1 + X2")

  parts           <- strsplit(f_char, "\\|")[[1]]
  lhs_rhs         <- trimws(parts[1])
  confounders_str <- trimws(parts[2])

  lhs_rhs_form <- stats::as.formula(lhs_rhs)
  mf  <- stats::model.frame(lhs_rhs_form, data = data, na.action = stats::na.omit)
  Y   <- stats::model.response(mf)
  Tr  <- mf[, 2]

  conf_form <- stats::as.formula(paste("~", confounders_str))
  X_mat     <- stats::model.matrix(conf_form,
                                   data = data[stats::complete.cases(data), ])[, -1, drop = FALSE]

  list(Y = as.numeric(Y), T = as.numeric(Tr), X = X_mat)
}


#' @keywords internal
.compute_ipw_weights <- function(Tr, gps_fit, family) {
  t_pred <- stats::fitted(gps_fit$model)
  t_res  <- stats::residuals(gps_fit$model)
  sigma  <- stats::sd(t_res)

  if (family == "gaussian") {
    f_t_given_x <- stats::dnorm(Tr, mean = t_pred, sd = sigma)
    f_t         <- stats::dnorm(Tr, mean = mean(Tr), sd = stats::sd(Tr))

  } else if (family == "gamma") {
    mu          <- t_pred
    phi         <- sigma^2 / mean(mu)
    shape       <- mu^2 / phi
    rate        <- mu / phi
    f_t_given_x <- stats::dgamma(Tr, shape = shape, rate = rate)
    f_t         <- stats::dgamma(Tr,
                                 shape = mean(Tr)^2 / stats::var(Tr),
                                 rate  = mean(Tr)  / stats::var(Tr))
  } else {
    Tr_sc       <- (Tr - min(Tr) + 1e-6) / (max(Tr) - min(Tr) + 2e-6)
    mu_sc       <- (t_pred - min(Tr) + 1e-6) / (max(Tr) - min(Tr) + 2e-6)
    phi         <- max(mu_sc * (1 - mu_sc) / sigma^2 - 1, 1)
    a           <- mu_sc * phi
    b           <- (1 - mu_sc) * phi
    f_t_given_x <- stats::dbeta(Tr_sc, a, b)
    a0          <- mean(Tr_sc) * phi
    b0          <- (1 - mean(Tr_sc)) * phi
    f_t         <- stats::dbeta(Tr_sc, a0, b0)
  }

  f_t / pmax(f_t_given_x, 1e-10)
}


#' @keywords internal
.fit_ipw_spline <- function(Y, T_basis, T_basis_grid, weights,
                            se_method, boot_reps, conf_level, Tr, t_grid) {
  df_wls <- data.frame(Y = Y, T_basis)
  fit    <- stats::lm(Y ~ ., data = df_wls, weights = weights)
  pred   <- as.numeric(cbind(1, T_basis_grid) %*% stats::coef(fit))

  if (se_method == "sandwich") {
    V    <- sandwich::vcovHC(fit, type = "HC3")
    X_g  <- cbind(1, T_basis_grid)
    se_g <- sqrt(pmax(rowSums((X_g %*% V) * X_g), 0))
  } else {
    boot_fn <- function(data, idx) {
      w_b  <- weights[idx]
      df_b <- data[idx, ]
      fb   <- stats::lm(Y ~ ., data = df_b, weights = w_b)
      as.numeric(cbind(1, T_basis_grid) %*% stats::coef(fb))
    }
    boot_df <- data.frame(Y = Y, T_basis)
    bt      <- boot::boot(boot_df, boot_fn, R = boot_reps)
    se_g    <- apply(bt$t, 2, stats::sd, na.rm = TRUE)
  }

  z     <- stats::qnorm(1 - (1 - conf_level) / 2)
  curve <- data.frame(
    t        = t_grid,
    estimate = pred,
    se       = se_g,
    lower    = pred - z * se_g,
    upper    = pred + z * se_g
  )
  list(curve = curve, fit = fit)
}


#' @keywords internal
.fit_gcomp_spline <- function(Y, Tr, T_basis, T_basis_grid, X,
                              spline_type, knots_exposure, t_range,
                              t_grid, se_method, boot_reps, conf_level) {
  df_ols <- data.frame(Y = Y, T_basis, X)
  fit    <- stats::lm(Y ~ ., data = df_ols)

  n         <- length(Y)
  ng        <- length(t_grid)
  pred_mean <- numeric(ng)
  spline_fn <- if (spline_type == "ns") splines::ns else splines::bs
  T_b_ref   <- spline_fn(Tr, knots = knots_exposure,
                         Boundary.knots = t_range, intercept = FALSE)

  for (j in seq_len(ng)) {
    T_b_j        <- stats::predict(T_b_ref, newx = rep(t_grid[j], n))
    df_j         <- data.frame(T_b_j, X)
    colnames(df_j) <- colnames(df_ols)[-1]
    pred_mean[j] <- mean(stats::predict(fit, newdata = df_j))
  }

  if (se_method == "bootstrap") {
    boot_fn <- function(data_list, idx) {
      Y_b  <- data_list$Y[idx]
      Tr_b <- data_list$Tr[idx]
      X_b  <- data_list$X[idx, , drop = FALSE]
      T_b  <- spline_fn(Tr_b, knots = knots_exposure,
                        Boundary.knots = t_range, intercept = FALSE)
      df_b <- data.frame(Y = Y_b, T_b, X_b)
      colnames(df_b) <- colnames(df_ols)
      fb   <- stats::lm(Y ~ ., data = df_b)
      sapply(t_grid, function(tv) {
        T_bj <- stats::predict(T_b, newx = rep(tv, length(Y_b)))
        dj   <- data.frame(T_bj, X_b)
        colnames(dj) <- colnames(df_ols)[-1]
        mean(stats::predict(fb, newdata = dj))
      })
    }
    data_list <- list(Y = Y, Tr = Tr, X = X)
    bt        <- boot::boot(data_list, boot_fn, R = boot_reps)
    se_g      <- apply(bt$t, 2, stats::sd, na.rm = TRUE)

  } else {
    V    <- sandwich::vcovHC(fit, type = "HC3")
    X_g  <- cbind(1, T_basis_grid,
                  matrix(colMeans(X), nrow = ng, ncol = ncol(X), byrow = TRUE))
    se_g <- sqrt(pmax(rowSums((X_g %*% V) * X_g), 0))
  }

  z     <- stats::qnorm(1 - (1 - conf_level) / 2)
  curve <- data.frame(
    t        = t_grid,
    estimate = pred_mean,
    se       = se_g,
    lower    = pred_mean - z * se_g,
    upper    = pred_mean + z * se_g
  )
  list(curve = curve, fit = fit)
}


#' @keywords internal
.combine_dr <- function(ipw, gcomp, conf_level) {
  est  <- 0.5 * ipw$curve$estimate + 0.5 * gcomp$curve$estimate
  se_g <- sqrt(0.5^2 * ipw$curve$se^2 + 0.5^2 * gcomp$curve$se^2)
  z    <- stats::qnorm(1 - (1 - conf_level) / 2)
  curve <- data.frame(
    t        = ipw$curve$t,
    estimate = est,
    se       = se_g,
    lower    = est - z * se_g,
    upper    = est + z * se_g
  )
  list(curve = curve)
}
