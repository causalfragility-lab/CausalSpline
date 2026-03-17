library(testthat)
library(CausalSpline)

# ── simulate_dose_response ───────────────────────────────────────────────────

test_that("simulate_dose_response returns correct structure", {
  dat <- simulate_dose_response(n = 200, dgp = "threshold", seed = 1)
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 200)
  expect_named(dat, c("Y", "T", "X1", "X2", "X3", "true_effect"))
  expect_true(all(dat$T >= 0 & dat$T <= 10))
})

test_that("all DGPs run without error", {
  dgps <- c("threshold", "diminishing", "nonmonotone", "linear", "sinusoidal")
  for (d in dgps) {
    expect_no_error(simulate_dose_response(100, dgp = d, seed = 1))
  }
})

# ── gps_model ────────────────────────────────────────────────────────────────

test_that("gps_model returns list with model, family, r_squared", {
  dat <- simulate_dose_response(300, seed = 2)
  X   <- as.matrix(dat[, c("X1", "X2", "X3")])
  res <- gps_model(dat$T, X, family = "gaussian")
  expect_type(res, "list")
  expect_named(res, c("model", "family", "r_squared"))
  expect_true(res$r_squared >= 0 && res$r_squared <= 1)
})

# ── trim_weights ─────────────────────────────────────────────────────────────

test_that("trim_weights clips extreme values", {
  w <- c(0.001, 1, 2, 3, 100)
  wt <- trim_weights(w, c(0.01, 0.99))
  expect_true(max(wt) < 100)
  expect_true(min(wt) > 0.001 || min(wt) == min(w))
})

# ── causal_spline: IPW ────────────────────────────────────────────────────────

test_that("causal_spline IPW returns correct class and curve", {
  dat <- simulate_dose_response(300, dgp = "linear", seed = 3)
  fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat,
                       method = "ipw", df_exposure = 3, eval_grid = 50)
  expect_s3_class(fit, "causal_spline")
  expect_equal(nrow(fit$curve), 50)
  expect_named(fit$curve, c("t", "estimate", "se", "lower", "upper"))
  expect_true(all(fit$curve$lower <= fit$curve$estimate))
  expect_true(all(fit$curve$upper >= fit$curve$estimate))
})

# ── causal_spline: G-computation ─────────────────────────────────────────────

test_that("causal_spline gcomp runs and produces sensible estimates", {
  dat <- simulate_dose_response(300, dgp = "diminishing", seed = 4)
  fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat,
                       method = "gcomp", df_exposure = 4, eval_grid = 40)
  expect_s3_class(fit, "causal_spline")
  expect_equal(nrow(fit$curve), 40)
  # Diminishing returns: curve should be broadly increasing over [0,10]
  coarse <- fit$curve[c(1, nrow(fit$curve)), "estimate"]
  expect_true(coarse[2] > coarse[1])
})

# ── predict.causal_spline ─────────────────────────────────────────────────────

test_that("predict.causal_spline returns data frame with correct rows", {
  dat  <- simulate_dose_response(300, seed = 5)
  fit  <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat,
                        method = "gcomp", df_exposure = 3)
  preds <- predict(fit, newt = c(1, 3, 5, 7, 9))
  expect_s3_class(preds, "data.frame")
  expect_equal(nrow(preds), 5)
  expect_named(preds, c("t", "estimate"))
})

# ── dose_response_curve ───────────────────────────────────────────────────────

test_that("dose_response_curve returns the curve data frame", {
  dat <- simulate_dose_response(200, seed = 6)
  fit <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat,
                       method = "ipw", eval_grid = 30)
  cr  <- dose_response_curve(fit)
  expect_s3_class(cr, "data.frame")
  expect_equal(nrow(cr), 30)
})

# ── Formula parsing error ─────────────────────────────────────────────────────

test_that("causal_spline errors with missing | in formula", {
  dat <- simulate_dose_response(100, seed = 7)
  expect_error(
    causal_spline(Y ~ T + X1 + X2, data = dat),
    "\\|"
  )
})
