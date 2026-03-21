#' Print method for causal_spline objects
#'
#' @param x A \code{causal_spline} object.
#' @param ... Ignored.
#'
#' @return Invisibly returns the input \code{causal_spline} object \code{x},
#'   unchanged. The function is called for its side effect of printing a
#'   compact summary to the console, showing the estimation method, spline
#'   type, degrees of freedom, sample size, treatment range, and number of
#'   evaluation points on the dose-response curve.
#'
#' @method print causal_spline
#' @export
print.causal_spline <- function(x, ...) {
  cat("CausalSpline fit\n")
  cat("  Method         :", toupper(x$method), "\n")
  cat("  Spline type    :", x$spline_type, "\n")
  cat("  df (exposure)  :", x$df_exposure, "\n")
  cat("  n              :", x$data_summary$n, "\n")
  cat("  Treatment range: [",
      round(x$t_range[1], 3), ",", round(x$t_range[2], 3), "]\n")
  cat("  Eval points    :", nrow(x$curve), "\n")
  invisible(x)
}


#' Summary method for causal_spline objects
#'
#' @param object A \code{causal_spline} object.
#' @param ... Ignored.
#'
#' @return Invisibly returns the input \code{causal_spline} object
#'   \code{object}, unchanged. The function is called for its side effect of
#'   printing a detailed summary to the console, including the original call,
#'   estimation method, spline configuration, treatment variable statistics,
#'   IPW diagnostics (effective sample size and weight range, if applicable),
#'   and a table of the estimated dose-response curve at seven representative
#'   percentile points (treatment value, point estimate, standard error, and
#'   confidence interval bounds).
#'
#' @method summary causal_spline
#' @export
summary.causal_spline <- function(object, ...) {
  cat("=== CausalSpline Summary ===\n\n")
  cat("Call:\n  "); print(object$call); cat("\n")
  cat("Estimation method  :", toupper(object$method), "\n")
  cat("Spline type        :", object$spline_type, "\n")
  cat("df (exposure)      :", object$df_exposure, "\n")
  cat("Interior knots     :", round(object$knots, 3), "\n\n")

  cat("Treatment variable:\n")
  cat("  n =", object$data_summary$n, "\n")
  cat("  Mean =", round(object$data_summary$t_mean, 3),
      "  SD =", round(object$data_summary$t_sd, 3), "\n")
  cat("  Range = [", round(object$t_range[1], 3),
      ",", round(object$t_range[2], 3), "]\n\n")

  if (!is.null(object$weights)) {
    ess <- sum(object$weights)^2 / sum(object$weights^2)
    cat("IPW diagnostics:\n")
    cat("  ESS =", round(ess), "/ n =", object$data_summary$n,
        "(", round(100 * ess / object$data_summary$n, 1), "%)\n")
    cat("  Weight range = [",
        round(min(object$weights), 3), ",",
        round(max(object$weights), 3), "]\n\n")
  }

  cat("Dose-response curve (selected percentiles):\n")
  idx <- round(seq(1, nrow(object$curve), length.out = 7))
  sub <- object$curve[idx, ]
  sub_print <- data.frame(
    t        = round(sub$t, 3),
    estimate = round(sub$estimate, 4),
    se       = round(sub$se, 4),
    lower    = round(sub$lower, 4),
    upper    = round(sub$upper, 4)
  )
  print(sub_print, row.names = FALSE)
  invisible(object)
}


#' Plot the estimated dose-response curve
#'
#' @description
#' Plots the estimated causal dose-response curve \eqn{E[Y(t)]} against \eqn{t}
#' with pointwise confidence bands and an optional rug for the observed
#' treatment distribution.
#'
#' @param x A \code{causal_spline} object.
#' @param rug Logical. Add a treatment distribution rug. Default \code{TRUE}.
#' @param truth Optional data frame with columns \code{t} and \code{true_effect}
#'   for overlaying the true dose-response (useful in simulations).
#' @param xlab Character. x-axis label. Default \code{"Treatment (T)"}.
#' @param ylab Character. y-axis label. Default \code{"E[Y(t)]"}.
#' @param title Character. Plot title. Default \code{NULL}.
#' @param ... Ignored.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' dat <- simulate_dose_response(200, dgp = "threshold")
#' fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat)
#' plot(fit)
#'
#' @method plot causal_spline
#' @export
plot.causal_spline <- function(x, rug = TRUE, truth = NULL,
                               xlab  = "Treatment (T)",
                               ylab  = "E[Y(t)]",
                               title = NULL, ...) {
  curve <- x$curve
  method_label <- switch(x$method,
                         ipw   = "IPW (GPS)",
                         gcomp = "G-computation",
                         dr    = "Doubly Robust"
  )
  if (is.null(title)) title <- paste("Causal Dose-Response Curve -", method_label)

  p <- ggplot2::ggplot(curve, ggplot2::aes(x = t, y = estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "#2166AC", alpha = 0.2) +
    ggplot2::geom_line(colour = "#2166AC", linewidth = 1.1) +
    ggplot2::geom_hline(yintercept = x$ate, linetype = "dashed",
                        colour = "grey50", linewidth = 0.7) +
    ggplot2::labs(title   = title, x = xlab, y = ylab,
                  caption = paste0(round(100 * x$conf_level),
                                   "% pointwise CI  |  ",
                                   "Dashed = marginal mean E[Y(t)]")) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  if (rug && !is.null(truth) && "T" %in% names(truth)) {
    p <- p + ggplot2::geom_rug(
      data = truth,
      ggplot2::aes(x = T, y = NULL),
      alpha = 0.25, sides = "b"
    )
  }

  if (!is.null(truth) && all(c("t", "true_effect") %in% names(truth))) {
    p <- p +
      ggplot2::geom_line(
        data = truth,
        ggplot2::aes(x = t, y = true_effect),
        colour = "#D6604D", linetype = "solid",
        linewidth = 1, inherit.aes = FALSE
      ) +
      ggplot2::labs(caption = paste0(
        round(100 * x$conf_level), "% pointwise CI  |  ",
        "Red = true curve  |  Dashed = E[Y(t)]"
      ))
  }
  p
}


#' Predict E[Y(t)] at new treatment values
#'
#' @param object A \code{causal_spline} object.
#' @param newt Numeric vector of treatment values at which to predict.
#' @param se_fit Logical. Return standard errors? Default \code{FALSE}.
#' @param warn_extrap Logical. Warn if any values in \code{newt} fall outside
#'   the observed treatment range? Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return A data frame with columns \code{t}, \code{estimate},
#'   \code{extrapolated}, and optionally \code{se}, \code{lower}, \code{upper}.
#'   The \code{extrapolated} column is \code{TRUE} for any \code{newt} value
#'   outside the observed treatment support.
#'
#' @examples
#' dat <- simulate_dose_response(200)
#' fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat)
#' predict(fit, newt = c(1, 2, 3, 4, 5))
#'
#' @method predict causal_spline
#' @export
predict.causal_spline <- function(object, newt, se_fit = FALSE,
                                  warn_extrap = TRUE, ...) {

  # ── Extrapolation check ────────────────────────────────────────────────────
  t_lo       <- object$t_range[1]
  t_hi       <- object$t_range[2]
  is_extrap  <- newt < t_lo | newt > t_hi

  if (warn_extrap && any(is_extrap)) {
    extrap_vals <- round(newt[is_extrap], 3)
    warning(
      "Some values in `newt` fall outside the observed treatment range [",
      round(t_lo, 3), ", ", round(t_hi, 3), "]: ",
      paste(extrap_vals, collapse = ", "), ".\n",
      "  Predictions at these points are extrapolations and may be unstable.",
      call. = FALSE
    )
  }

  # ── Spline basis at new points ─────────────────────────────────────────────
  spline_fn   <- if (object$spline_type == "ns") splines::ns else splines::bs
  T_basis_new <- stats::predict(
    spline_fn(seq(t_lo, t_hi, length.out = 100),
              knots          = object$knots,
              Boundary.knots = object$t_range,
              intercept      = FALSE),
    newx = newt
  )

  # ── Build design matrix ────────────────────────────────────────────────────
  coefs <- stats::coef(object$outcome_fit)
  n_sp  <- ncol(T_basis_new)
  n_cov <- length(coefs) - 1 - n_sp
  if (n_cov > 0) {
    X_mean <- matrix(0, nrow = length(newt), ncol = n_cov)
    Xmat   <- cbind(1, T_basis_new, X_mean)
  } else {
    Xmat   <- cbind(1, T_basis_new)
  }

  est <- as.numeric(Xmat %*% coefs[seq_len(ncol(Xmat))])
  out <- data.frame(t = newt, estimate = est, extrapolated = is_extrap)

  # ── Optional SEs ──────────────────────────────────────────────────────────
  if (se_fit) {
    V    <- sandwich::vcovHC(object$outcome_fit, type = "HC3")
    se_v <- sqrt(pmax(rowSums((Xmat %*% V[seq_len(ncol(Xmat)),
                                          seq_len(ncol(Xmat))]) * Xmat), 0))
    z         <- stats::qnorm(1 - (1 - object$conf_level) / 2)
    out$se    <- se_v
    out$lower <- est - z * se_v
    out$upper <- est + z * se_v
  }

  out
}


#' Extract the dose-response curve data frame
#'
#' @description
#' Convenience function to pull the estimated \eqn{E[Y(t)]} curve from a
#' fitted \code{causal_spline} object.
#'
#' @param fit A \code{causal_spline} object.
#'
#' @return A data frame with columns \code{t}, \code{estimate}, \code{se},
#'   \code{lower}, \code{upper}.
#'
#' @examples
#' dat <- simulate_dose_response(200)
#' fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat)
#' curve_df <- dose_response_curve(fit)
#' head(curve_df)
#'
#' @export
dose_response_curve <- function(fit) {
  stopifnot(inherits(fit, "causal_spline"))
  fit$curve
}
