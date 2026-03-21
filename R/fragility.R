#' Numerical derivatives of a CausalSpline dose-response curve
#'
#' @description
#' Computes first and second numerical derivatives of the fitted dose-response
#' curve \eqn{E[Y(t)]} using central finite differences applied to
#' \code{\link{predict.causal_spline}}. Useful for identifying regions of rapid
#' change (high first derivative) or inflection / diminishing returns (second
#' derivative changes sign).
#'
#' @param object A fitted \code{causal_spline} object.
#' @param grid Optional numeric vector of treatment values at which to evaluate
#'   the derivatives. If \code{NULL}, the fitted evaluation grid stored in
#'   \code{object$curve$t} is used.
#' @param h Numeric. Step size for finite differences. If \code{NULL}, chosen
#'   automatically as \eqn{(t_{max} - t_{min}) / 500}.
#' @param eps Small positive constant used by \code{\link{fragility_curve}} to
#'   stabilise division. Default \code{1e-6}.
#' @param ... Ignored.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{t}}{Treatment grid values.}
#'   \item{\code{estimate}}{Fitted \eqn{E[Y(t)]}.}
#'   \item{\code{se}}{Standard error of \eqn{E[Y(t)]} from the fitted curve.}
#'   \item{\code{derivative}}{First derivative \eqn{f'(t)}.}
#'   \item{\code{second_derivative}}{Second derivative \eqn{f''(t)}.}
#' }
#'
#' @examples
#' dat <- simulate_dose_response(200, dgp = "threshold")
#' fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat)
#' gd  <- gradient_curve(fit)
#' head(gd)
#'
#' @export
gradient_curve <- function(object, grid = NULL, h = NULL, eps = 1e-6, ...) {
  if (!inherits(object, "causal_spline"))
    stop("`object` must be a fitted causal_spline object.")

  # ── Evaluation grid ────────────────────────────────────────────────────────
  if (is.null(grid)) grid <- object$curve$t
  grid <- sort(unique(as.numeric(grid)))

  # ── Treatment range ────────────────────────────────────────────────────────
  t_lo <- object$t_range[1]
  t_hi <- object$t_range[2]

  # ── Step size ──────────────────────────────────────────────────────────────
  if (is.null(h)) h <- (t_hi - t_lo) / 500
  if (h <= 0) stop("`h` must be positive.")

  # ── Clamped left / right evaluation points ─────────────────────────────────
  left  <- pmax(grid - h, t_lo)
  right <- pmin(grid + h, t_hi)

  # ── Fitted values + SE (suppress extrapolation warnings) ──────────────────
  pred0 <- predict(object, newt = grid,  se_fit = TRUE,  warn_extrap = FALSE)
  predl <- predict(object, newt = left,  se_fit = FALSE, warn_extrap = FALSE)
  predr <- predict(object, newt = right, se_fit = FALSE, warn_extrap = FALSE)

  mu0 <- pred0$estimate
  se0 <- pred0$se
  mul <- predl$estimate
  mur <- predr$estimate

  # ── Effective step sizes after clamping ────────────────────────────────────
  hl <- grid - left
  hr <- right - grid

  # ── First derivative ───────────────────────────────────────────────────────
  d1         <- rep(NA_real_, length(grid))
  interior   <- hl > 0 & hr > 0
  left_only  <- hl == 0 & hr > 0
  right_only <- hl > 0 & hr == 0

  d1[interior]   <- (mur[interior] - mul[interior]) /
    (hr[interior] + hl[interior])
  d1[left_only]  <- (mur[left_only] - mu0[left_only]) / hr[left_only]
  d1[right_only] <- (mu0[right_only] - mul[right_only]) / hl[right_only]

  # ── Second derivative ──────────────────────────────────────────────────────
  d2        <- rep(NA_real_, length(grid))
  good2     <- hl > 0 & hr > 0
  d2[good2] <- 2 * (
    mul[good2] / (hl[good2] * (hl[good2] + hr[good2])) -
      mu0[good2] / (hl[good2] * hr[good2]) +
      mur[good2] / (hr[good2] * (hl[good2] + hr[good2]))
  )

  data.frame(
    t                 = grid,
    estimate          = mu0,
    se                = se0,
    derivative        = d1,
    second_derivative = d2
  )
}


#' Geometric fragility curve for a CausalSpline fit
#'
#' @description
#' Computes a pointwise shape-based fragility measure from the first and second
#' derivatives of the fitted dose-response curve, with five enhancements:
#'
#' \enumerate{
#'   \item \strong{Adaptive eps}: stabilisation constant is
#'     \eqn{0.05 \times \text{median}(|f'(t)|)} so interpretation is
#'     consistent across datasets.
#'   \item \strong{Interpretation zones}: fragility is classified into
#'     \code{"low"}, \code{"moderate"}, and \code{"high"} based on the 50th
#'     and 75th percentiles of the pointwise fragility distribution.
#'   \item \strong{Uncertainty-normalised fragility}: an additional column
#'     \eqn{F^*(t) = F(t) / \text{SE}(\hat{\mu}(t))} combines shape
#'     instability with statistical uncertainty.
#'   \item \strong{Support density}: the kernel density of the treatment
#'     variable (if \code{t_obs} is supplied) is attached, flagging regions
#'     with sparse data.
#'   \item \strong{High-fragility flag}: logical column \code{high_fragility}
#'     marks points above the 75th percentile.
#' }
#'
#' @param object A fitted \code{causal_spline} object.
#' @param grid Optional numeric evaluation grid. If \code{NULL}, the fitted
#'   grid in \code{object$curve$t} is used.
#' @param h Step size for finite differences. Default \code{NULL} (auto).
#' @param eps Numeric or \code{NULL}. If \code{NULL} (default), set adaptively
#'   to \eqn{0.05 \times \text{median}(|f'(t)|)}. If numeric, used directly.
#' @param type Character. Fragility definition: \code{"curvature_ratio"}
#'   (default) or \code{"inverse_slope"}.
#' @param t_obs Optional numeric vector of observed treatment values (used to
#'   compute support density). If \code{NULL}, density column is omitted.
#' @param ... Ignored.
#'
#' @return A data frame of class \code{"fragility_curve"} with columns:
#' \describe{
#'   \item{\code{t}}{Treatment grid.}
#'   \item{\code{estimate}}{Fitted \eqn{E[Y(t)]}.}
#'   \item{\code{se}}{Standard error of fitted curve.}
#'   \item{\code{derivative}}{First derivative.}
#'   \item{\code{second_derivative}}{Second derivative.}
#'   \item{\code{fragility}}{Pointwise fragility \eqn{F(t)}.}
#'   \item{\code{fragility_norm}}{Uncertainty-normalised fragility
#'     \eqn{F^*(t) = F(t) / \text{SE}(\hat{\mu}(t))}.}
#'   \item{\code{fragility_zone}}{Factor: \code{"low"}, \code{"moderate"},
#'     \code{"high"}.}
#'   \item{\code{high_fragility}}{Logical: \code{TRUE} if above 75th
#'     percentile.}
#'   \item{\code{support_density}}{Kernel density of T at each grid point
#'     (only if \code{t_obs} supplied).}
#'   \item{\code{fragility_type}}{Character. Type used.}
#' }
#'
#' @examples
#' dat <- simulate_dose_response(200, dgp = "threshold")
#' fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat)
#' fc  <- fragility_curve(fit, t_obs = dat$T)
#' plot(fc)
#'
#' @export
fragility_curve <- function(object,
                            grid  = NULL,
                            h     = NULL,
                            eps   = NULL,
                            type  = c("curvature_ratio", "inverse_slope"),
                            t_obs = NULL,
                            ...) {
  type <- match.arg(type)

  g <- gradient_curve(object = object, grid = grid, h = h, ...)

  # ── (A) Adaptive eps ───────────────────────────────────────────────────────
  if (is.null(eps)) {
    med_slope <- stats::median(abs(g$derivative), na.rm = TRUE)
    eps       <- 0.05 * max(med_slope, 1e-8)
  }

  # ── Raw fragility ──────────────────────────────────────────────────────────
  g$fragility <- if (type == "curvature_ratio") {
    abs(g$second_derivative) / (abs(g$derivative) + eps)
  } else {
    1 / (abs(g$derivative) + eps)
  }

  # ── (E) Uncertainty-normalised fragility F*(t) = F(t) / SE(mu(t)) ─────────
  g$fragility_norm <- g$fragility / (g$se + 1e-10)

  # ── (B) Interpretation zones ───────────────────────────────────────────────
  q50 <- stats::quantile(g$fragility, 0.50, na.rm = TRUE)
  q75 <- stats::quantile(g$fragility, 0.75, na.rm = TRUE)

  g$high_fragility <- !is.na(g$fragility) & g$fragility > q75
  g$fragility_zone <- ifelse(
    is.na(g$fragility), NA_character_,
    ifelse(g$fragility > q75, "high",
           ifelse(g$fragility > q50, "moderate", "low"))
  )
  g$fragility_zone <- factor(g$fragility_zone,
                             levels = c("low", "moderate", "high"))

  # ── (D) Support density ────────────────────────────────────────────────────
  if (!is.null(t_obs)) {
    dens              <- stats::density(t_obs, from = min(g$t),
                                        to   = max(g$t),
                                        n    = length(g$t))
    g$support_density <- stats::approx(dens$x, dens$y, xout = g$t)$y
  }

  g$fragility_type <- type
  class(g) <- c("fragility_curve", "data.frame")
  g
}


#' Regional fragility summary for a CausalSpline fit
#'
#' @description
#' Integrates the pointwise fragility curve over a treatment interval
#' \eqn{[a, b]} using the trapezoidal rule. Useful for comparing sensitivity
#' across dose ranges (e.g., low vs. high dose) or summarising instability
#' at a policy-relevant threshold.
#'
#' @param object A fitted \code{causal_spline} object.
#' @param a Numeric scalar. Lower bound of the integration interval.
#' @param b Numeric scalar. Upper bound of the integration interval.
#' @param grid Optional numeric evaluation grid within \eqn{[a, b]}. If
#'   \code{NULL}, a dense grid of 300 points is created automatically.
#' @param h Step size for finite differences. Default \code{NULL} (auto).
#' @param eps Adaptive eps passed to \code{\link{fragility_curve}}.
#'   Default \code{NULL} (auto).
#' @param type Fragility definition: \code{"curvature_ratio"} (default) or
#'   \code{"inverse_slope"}.
#' @param normalize Logical. Divide integral by interval length? Default
#'   \code{TRUE}.
#' @param t_obs Optional numeric vector of observed treatment values for
#'   support density. Default \code{NULL}.
#' @param ... Ignored.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{interval}}{Named numeric vector \code{c(a, b)} after clamping
#'     to the observed support.}
#'   \item{\code{type}}{Fragility type used.}
#'   \item{\code{integral_fragility}}{Trapezoidal integral of fragility over
#'     \eqn{[a, b]}.}
#'   \item{\code{average_fragility}}{Integral divided by interval length.}
#'   \item{\code{normalized}}{Logical flag.}
#'   \item{\code{curve}}{The full \code{fragility_curve} data frame.}
#' }
#'
#' @examples
#' dat <- simulate_dose_response(200, dgp = "threshold")
#' fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat)
#' region_fragility(fit, a = 2, b = 5)
#'
#' @export
region_fragility <- function(object,
                             a,
                             b,
                             grid      = NULL,
                             h         = NULL,
                             eps       = NULL,
                             type      = c("curvature_ratio", "inverse_slope"),
                             normalize = TRUE,
                             t_obs     = NULL,
                             ...) {
  type <- match.arg(type)

  if (!is.numeric(a) || !is.numeric(b) ||
      length(a) != 1  || length(b) != 1)
    stop("`a` and `b` must be numeric scalars.")
  if (a >= b)
    stop("Require a < b.")

  t_lo <- object$t_range[1]
  t_hi <- object$t_range[2]
  a0   <- max(a, t_lo)
  b0   <- min(b, t_hi)

  if (a0 >= b0)
    stop("Requested interval does not overlap the observed treatment support.")

  if (is.null(grid)) {
    grid <- seq(a0, b0, length.out = 300)
  } else {
    grid <- sort(unique(as.numeric(grid)))
    grid <- grid[grid >= a0 & grid <= b0]
    if (length(grid) < 2)
      stop("Grid within [a, b] must contain at least 2 points.")
  }

  fc <- fragility_curve(object = object, grid = grid, h = h,
                        eps = eps, type = type, t_obs = t_obs, ...)

  ok <- is.finite(fc$t) & is.finite(fc$fragility)
  x  <- fc$t[ok]
  y  <- fc$fragility[ok]

  if (length(x) < 2)
    stop("Not enough finite points to integrate fragility over the interval.")

  integral      <- sum(diff(x) * (utils::head(y, -1) +
                                    utils::tail(y, -1)) / 2)
  avg_fragility <- integral / (max(x) - min(x))

  list(
    interval           = c(a = a0, b = b0),
    type               = type,
    integral_fragility = integral,
    average_fragility  = avg_fragility,
    normalized         = normalize,
    curve              = fc
  )
}


#' Plot method for fragility_curve objects
#'
#' @description
#' Produces a dual-panel plot: the dose-response curve (top) and the
#' fragility curve (bottom), with high-fragility regions shaded. If
#' \code{support_density} is present in \code{x} (i.e. \code{t_obs} was
#' supplied to \code{\link{fragility_curve}}), a scaled density ribbon is
#' overlaid on the fragility panel to flag low-support regions.
#'
#' @param x A \code{fragility_curve} data frame (output of
#'   \code{\link{fragility_curve}}).
#' @param ... Ignored.
#'
#' @return A combined \code{patchwork} plot if the \pkg{patchwork} package is
#'   installed, otherwise the two panels are printed separately and a list of
#'   two \code{ggplot2} objects is returned invisibly.
#'
#' @examples
#' dat <- simulate_dose_response(200, dgp = "threshold")
#' fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat)
#' fc  <- fragility_curve(fit, t_obs = dat$T)
#' plot(fc)
#'
#' @method plot fragility_curve
#' @export
plot.fragility_curve <- function(x, ...) {

  ftype <- unique(x$fragility_type)
  q75   <- stats::quantile(x$fragility, 0.75, na.rm = TRUE)
  q50   <- stats::quantile(x$fragility, 0.50, na.rm = TRUE)

  # ── Top panel: dose-response curve with fragility shading ─────────────────
  hi <- x[!is.na(x$high_fragility) & x$high_fragility, ]

  p_curve <- ggplot2::ggplot(x, ggplot2::aes(x = t, y = estimate)) +
    {
      if (nrow(hi) > 0)
        ggplot2::annotate("rect",
                          xmin = min(hi$t), xmax = max(hi$t),
                          ymin = -Inf,      ymax = Inf,
                          fill = "#D6604D", alpha = 0.08)
      else
        NULL
    } +
    ggplot2::geom_line(colour = "#2166AC", linewidth = 1.1) +
    ggplot2::labs(
      title    = "Dose-response curve with fragility regions",
      subtitle = paste0("Shaded = high fragility (top 25%)  |  Type: ", ftype),
      x        = NULL,
      y        = "E[Y(t)]"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold"),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  # ── Bottom panel: fragility curve ─────────────────────────────────────────
  p_frag <- ggplot2::ggplot(x, ggplot2::aes(x = t, y = fragility)) +
    ggplot2::geom_hline(yintercept = q75, linetype = "dashed",
                        colour = "#D6604D", linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = q50, linetype = "dotted",
                        colour = "grey60", linewidth = 0.6) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ifelse(fragility > q75, q75,       NA_real_),
                   ymax = ifelse(fragility > q75, fragility, NA_real_)),
      fill = "#D6604D", alpha = 0.25, na.rm = TRUE
    ) +
    ggplot2::geom_line(colour = "#2166AC", linewidth = 0.9, na.rm = TRUE) +
    ggplot2::labs(
      x       = "Treatment (T)",
      y       = "Fragility",
      caption = paste0("Red dashed = 75th pct (high)  |  ",
                       "Grey dotted = 50th pct (moderate)")
    ) +
    ggplot2::theme_bw(base_size = 12)

  # ── (D) Support density overlay ───────────────────────────────────────────
  if ("support_density" %in% names(x)) {
    frag_range    <- range(x$fragility, na.rm = TRUE)
    dens_max      <- max(x$support_density, na.rm = TRUE)
    scale_fac     <- diff(frag_range) * 0.3 / max(dens_max, 1e-10)
    # Compute scaled density and attach to a local copy of x
    x_dens        <- x
    x_dens$dens_scaled <- x$support_density * scale_fac + frag_range[1]

    p_frag <- p_frag +
      ggplot2::geom_area(
        data         = x_dens,
        ggplot2::aes(x = t, y = dens_scaled),
        fill         = "grey70",
        alpha        = 0.35,
        na.rm        = TRUE,
        inherit.aes  = FALSE
      ) +
      ggplot2::labs(caption = paste0(
        "Red dashed = 75th pct  |  ",
        "Grey area = treatment support density"
      ))
  }

  # ── Combine panels ─────────────────────────────────────────────────────────
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p_curve, p_frag, ncol = 1)
  } else {
    print(p_curve)
    print(p_frag)
    invisible(list(curve = p_curve, fragility = p_frag))
  }
}
