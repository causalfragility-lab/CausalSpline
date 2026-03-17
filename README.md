# CausalSpline <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/CausalSpline)](https://CRAN.R-project.org/package=CausalSpline)
[![R-CMD-check](https://github.com/yourgithub/CausalSpline/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourgithub/CausalSpline/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

**CausalSpline** estimates **nonlinear causal dose-response functions** for continuous treatments using spline-based methods under standard causal assumptions (unconfoundedness + positivity).

Most causal inference tools assume a linear treatment effect:

$$Y = \beta_0 + \beta_1 T + \gamma X + \varepsilon$$

Real policy and health problems often violate this: dosage effects, thresholds, diminishing returns, and non-monotone relationships are common. CausalSpline models the exposure-response curve nonparametrically:

$$E[Y(t)] = \beta_0 + f(T)$$

where $f(T)$ is a natural cubic spline or B-spline, estimated via **IPW**, **G-computation**, or **doubly-robust** methods.

---

## Installation

```r
# CRAN (once published)
install.packages("CausalSpline")

# Development version
remotes::install_github("yourgithub/CausalSpline")
```

---

## Quick example

```r
library(CausalSpline)

# Simulate data with a threshold effect
set.seed(42)
dat <- simulate_dose_response(n = 600, dgp = "threshold", confounding = 0.6)

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

# Predict at specific treatment values
predict(fit_gc, newt = c(1, 3, 5, 7, 9))

# Check positivity / overlap
check_overlap(dat$T, fit$weights)
```

---

## Supported methods

| Method | `method =` | Consistent if ... |
|--------|------------|-------------------|
| Inverse Probability Weighting | `"ipw"` | GPS model correct |
| G-computation | `"gcomp"` | Outcome model correct |
| Doubly Robust (AIPW) | `"dr"` | At least one model correct |

---

## Supported DGPs (for simulation)

```r
simulate_dose_response(dgp = "threshold")    # flat, then linear rise
simulate_dose_response(dgp = "diminishing")  # concave / log shape
simulate_dose_response(dgp = "nonmonotone")  # inverted-U
simulate_dose_response(dgp = "linear")       # linear (baseline)
simulate_dose_response(dgp = "sinusoidal")   # oscillatory (hard)
```

---

## Relationship to MultiSpline

CausalSpline is the **causal extension** of [MultiSpline](https://CRAN.R-project.org/package=MultiSpline). While MultiSpline fits multivariate spline regressions for prediction, CausalSpline:

- Incorporates the **generalised propensity score** for causal identification
- Targets the **marginal structural mean** $E[Y(t)]$, not the conditional mean
- Provides **doubly-robust** estimation
- Includes **positivity diagnostics** (ESS, weight trimming)

---

## References

- Hirano, K., & Imbens, G. W. (2004). *The propensity score with continuous treatments.*  
- Imbens, G. W. (2000). *The role of the propensity score in estimating dose-response functions.* Biometrika.  
- Flores et al. (2012). *Estimating the effects of length of exposure to instruction.* Rev. Econ. Stat.

---

## Citation

```r
citation("CausalSpline")
```
