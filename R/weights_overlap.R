#' Trim extreme IPW weights
#'
#' @description
#' Winsorises (clips) weights at specified quantiles to reduce variance from
#' extreme propensity scores. This is a standard stabilisation step when
#' using inverse probability weighting for continuous treatments.
#'
#' @param weights Numeric vector of (possibly unstabilised) IPW weights.
#' @param quantiles Numeric vector of length 2: lower and upper quantile
#'   thresholds. Default \code{c(0.01, 0.99)}.
#'
#' @return Numeric vector of trimmed weights (same length as input).
#'
#' @examples
#' w <- rexp(200)
#' w_trimmed <- trim_weights(w, c(0.01, 0.99))
#' summary(w_trimmed)
#'
#' @export
trim_weights <- function(weights, quantiles = c(0.01, 0.99)) {
  stopifnot(length(quantiles) == 2, quantiles[1] < quantiles[2])
  lo <- stats::quantile(weights, quantiles[1], na.rm = TRUE)
  hi <- stats::quantile(weights, quantiles[2], na.rm = TRUE)
  pmin(pmax(weights, lo), hi)
}


#' Diagnose positivity / overlap for a continuous treatment
#'
#' @description
#' Plots the distribution of the treatment variable conditional on covariate
#' strata and returns effective sample size (ESS) and weight diagnostics to
#' assess the positivity (overlap) assumption.
#'
#' @param treatment Numeric vector of treatment values.
#' @param weights Numeric vector of IPW weights (length must equal
#'   \code{length(treatment)}).
#' @param plot Logical. If \code{TRUE}, returns a \code{ggplot2} diagnostic
#'   plot. Default \code{TRUE}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{ess}}{Effective sample size: \eqn{(\sum w_i)^2 / \sum w_i^2}.}
#'   \item{\code{weight_summary}}{Five-number summary of the weights.}
#'   \item{\code{plot}}{ggplot2 object (if \code{plot = TRUE}).}
#' }
#'
#' @examples
#' dat <- simulate_dose_response(200)
#' X   <- as.matrix(dat[, c("X1", "X2", "X3")])
#' gps <- gps_model(dat$T, X)
#' w   <- trim_weights(abs(1 / stats::residuals(gps$model)), c(0.01, 0.99))
#' check_overlap(dat$T, w)
#'
#' @export
check_overlap <- function(treatment, weights, plot = TRUE) {
  n   <- length(treatment)
  ess <- sum(weights)^2 / sum(weights^2)
  ws  <- summary(weights)

  p <- NULL
  if (plot) {
    df_plot <- data.frame(t = treatment, w = weights)
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = t, weight = w)) +
      ggplot2::geom_histogram(bins = 30, fill = "#2166AC", alpha = 0.7,
                              colour = "white") +
      ggplot2::geom_rug(ggplot2::aes(x = t, y = NULL),
                        alpha = 0.3, sides = "b") +
      ggplot2::labs(
        title    = "Weighted treatment distribution (overlap diagnostic)",
        subtitle = paste0("ESS = ", round(ess), " / n = ", n,
                          " (", round(100 * ess / n, 1), "%)"),
        x = "Treatment T",
        y = "Weighted count"
      ) +
      ggplot2::theme_bw()
  }

  list(ess = ess, weight_summary = ws, plot = p)
}
