# JoCuR

> **Joint Cure Rate modeling for longitudinal tumor burden and time-to-event outcomes**

[![R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status](https://img.shields.io/badge/status-research%20ready-brightgreen)]()

`JoCuR` implements a C++-accelerated joint cure rate framework for oncology-style data, where tumor burden trajectories and progression times are modeled together under latent subgroup structure.

## Why JoCuR?

Many oncology longitudinal profiles are **not purely linear**:

- some patients show **change-point dynamics** (decrease then regrowth),
- others remain **stable/long-term controlled**,
- censoring and event timing carry critical information about latent trajectory shape.

`JoCuR` is designed exactly for this setting.

## Key Features

- ⚡ **C++-accelerated MCEM** estimation for practical runtime.
- 🧠 **Joint modeling** of longitudinal outcome + time-to-event.
- 🎯 **Cure/stable subgroup modeling** with latent membership.
- 🧪 **Built-in simulation** for benchmarking and method studies.
- 🔁 **Two inference paths** via `non_std`:
  - `non_std = 0` (default): standardized identifiable implementation,
  - `non_std = 1`: legacy non-standard path.

## Installation

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("quyixiang/JoCuR")
```

## Quick Start

```r
library(JoCuR)

set.seed(123)
fmla.tte  <- as.formula(survival::Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE)
fmla.long <- as.formula(PCHG ~ 0 + Y0SCALE)

sim_data <- CRsimulation(
  fmla.tte = fmla.tte,
  fmla.long = fmla.long,
  beta.tte = c(0.2), normal.tte = TRUE, sd.tte = 0.2,
  beta.y = c(-0.5), sd.y = 0.1,
  beta.y.cure = c(-0.2), sd.y.cured = 0.2,
  cured.rate = 0.4,
  n = 100,
  censor.parameter = 0.5,
  time.interval = 0.1
)

survdat <- sim_data[["survdat"]]
longdat <- sim_data[["longdat"]]

init <- list(
  mu_r = c(0.5, 0, -0.5, 0.5),
  Sigma_r = diag(rep(0.2, 4)),
  beta_y = c(-0.5), sigma_y_sq = 0.01,
  beta_tte = c(0.2), sigma_tte_sq = 0.04,
  pi_c = 0.4,
  beta_cure = c(-0.2),
  mu_r_cure = c(0, -0.2),
  Sigma_r_cure = diag(c(0.2, 0.2)),
  sigma_y_cure_sq = 0.04
)

fit <- infer_CRJoint_MLE(
  survdat = survdat,
  longdat = longdat,
  init = init,
  fmla.tte = fmla.tte,
  fmla.long = fmla.long,
  longdat.time = "visittime",
  id.indicator = "id",
  maxIter = 50,
  tol = 5e-3,
  non_std = 0
)

ey_t05 <- E_y(
  X = survdat$Y0SCALE,
  W = survdat$Y0SCALE,
  cure_res = fit,
  t_new = 0.5,
  J = 100
)

ey_t05
```

## Main API

- `CRsimulation(...)`  
  Generate simulation datasets (`survdat`, `longdat`) with cure/stable structure.

- `infer_CRJoint_MLE(...)`  
  Fit the joint cure rate model and return estimated parameters.

- `E_y(X, W, cure_res, t_new, J)`  
  Estimate expected longitudinal outcome at a chosen time point.

## Project Status

- Package checks pass locally via `R CMD build` and `R CMD check`.
- Suitable for research workflows and method development.

## Citation

If this package supports your work, please cite:

Qu Y, et al. *Cure Rate Joint Model for Time-to-Event Data and Longitudinal Tumor Burden with Potential Change Points*.

## Author

**Yixiang Qu**  
Email: `yqu@unc.edu`
