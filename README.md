# CausalSpline <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/CausalSpline)](https://CRAN.R-project.org/package=CausalSpline)
[![R-CMD-check](https://github.com/causalfragility-lab/CausalSpline/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/causalfragility-lab/CausalSpline/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

**CausalSpline** estimates **nonlinear causal dose-response functions** for
continuous treatments using spline-based methods under standard causal
assumptions (unconfoundedness + positivity).

Most causal inference tools assume a linear treatment effect:

$$Y = \beta_0 + \beta_1 T + \gamma X + \varepsilon$$

Real policy and health problems often violate this: dosage effects, thresholds,
diminishing returns, and non-monotone relationships are common. CausalSpline
models the exposure-response curve nonparametrically:

$$E[Y(t)] = \beta_0 + f(T)$$

where $f(T)$ is a natural cubic spline or B-spline, estimated via **IPW**,
**G-computation**, or **doubly-robust** methods.

---

## Installation

```r
# CRAN
install.packages("CausalSpline")

# Development version
remotes::install_github("causalfragility-lab/CausalSpline")
```

---

## Quick Example

```r
library(CausalSpline)

# Simulate data with a threshold effect
set.seed(42)
dat <- simulate_dose_response(n = 500, dgp = "threshold", confounding = 0.6)

# Fit causal dose-response curve via IPW
fit <- causal_spline(
  Y ~ T | X1 + X2 + X3,   # outcome ~ treatment | confounders
  data        = dat,
  method      = "ipw",
  df_exposure = 5
)

summary(fit)
plot(fit)

# G-computation alternative
fit_gc <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat,
                        method = "gcomp", df_exposure = 5)

# Doubly robust estimation
fit_dr <- causal_spline(Y ~ T | X1 + X2 + X3, data = dat,
                        method = "dr", df_exposure = 5)

# Predict at specific treatment values
predict(fit_gc, newt = c(3, 5, 7), se_fit = TRUE)

# Check positivity / overlap
check_overlap(dat$T, fit$weights)
```

---

## Fragility Diagnostics

A key feature of CausalSpline is **geometric fragility diagnostics** — tools
for evaluating the structural stability of the estimated dose-response curve
across treatment levels.

```r
# Slope-based fragility (identifies flat / weak-effect regions)
fc_slope <- fragility_curve(fit, type = "inverse_slope")
plot(fc_slope)

# Curvature-based fragility (identifies thresholds and turning points)
fc_curv <- fragility_curve(fit, type = "curvature_ratio")
plot(fc_curv)

# Regional fragility summary over a policy-relevant interval
region_fragility(fit, a = 3, b = 6, type = "curvature_ratio")
```

These diagnostics complement traditional sensitivity analyses by examining the
**internal geometric stability** of the causal function — where is the estimated
relationship flat, rapidly changing, or structurally ambiguous?

---

## Supported Methods

| Method | `method =` | Consistent if ... |
|---|---|---|
| Inverse Probability Weighting | `"ipw"` | GPS model correct |
| G-computation | `"gcomp"` | Outcome model correct |
| Doubly Robust (AIPW) | `"dr"` | At least one model correct |

---

## Supported DGPs (for Simulation)

```r
simulate_dose_response(dgp = "threshold")   # flat region then linear rise
simulate_dose_response(dgp = "diminishing") # concave / diminishing returns
simulate_dose_response(dgp = "nonmonotone") # inverted-U / hump shape
simulate_dose_response(dgp = "linear")      # linear baseline
simulate_dose_response(dgp = "sinusoidal")  # oscillatory / complex shape
```

---

## Key Functions

| Function | Description |
|---|---|
| `causal_spline()` | Estimate the causal dose-response curve |
| `fragility_curve()` | Compute local geometric fragility diagnostics |
| `region_fragility()` | Regional fragility summary over treatment intervals |
| `gradient_curve()` | First and second derivatives of the estimated curve |
| `dose_response_curve()` | Extract the dose-response data frame |
| `predict()` | Predict E[Y(t)] at new treatment values |
| `check_overlap()` | Overlap / positivity diagnostics (ESS, weight plot) |
| `simulate_dose_response()` | Simulate nonlinear dose-response datasets |

---

## Relationship to MultiSpline

CausalSpline is the **causal extension** of
[MultiSpline](https://CRAN.R-project.org/package=MultiSpline). While
MultiSpline fits multivariate spline regressions for prediction, CausalSpline:

- Incorporates the **generalised propensity score** for causal identification
- Targets the **marginal structural mean** $E[Y(t)]$, not the conditional mean
- Provides **doubly-robust** estimation
- Includes **positivity diagnostics** (ESS, weight trimming)
- Adds **geometric fragility diagnostics** with no precedent in existing packages

---

## References

- Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous
  treatments. In *Applied Bayesian Modeling and Causal Inference from
  Incomplete-Data Perspectives* (pp. 73–84). Wiley.
- Imbens, G. W. (2000). The role of the propensity score in estimating
  dose-response functions. *Biometrika*, 87(3), 706–710.
- Robins, J. M., Hernán, M. A., & Brumback, B. (2000). Marginal structural
  models and causal inference in epidemiology. *Epidemiology*, 11(5), 550–560.
- Flores, C. A., Flores-Lagunes, A., Gonzalez, A., & Neumann, T. C. (2012).
  Estimating the effects of length of exposure to instruction in a training
  program. *Review of Economics and Statistics*, 94(1), 153–171.

---

## Citation

```r
citation("CausalSpline")
```

---

## License

GPL (>= 3) — see [LICENSE.md](LICENSE.md) for details.
 
Maintainer: Subir Hait <haitsubi@msu.edu>  
Michigan State University
