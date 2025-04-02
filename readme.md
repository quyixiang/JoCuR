# JoCuR: Joint Cure Rate Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**JoCuR** is an R package that implements a C++-accelerated joint cure rate model for longitudinal tumor burden and time-to-event data, as described in the paper *"Change-Point Detection Using a Cure Rate Joint Model for Longitudinal Tumor Burden and Time-to-Event Data"* by Yixiang Qu.

- **Author**: Yixiang Qu (<yqu@unc.edu>)
- **License**: MIT

## Description

The `JoCuR` package implements a C++-accelerated joint cure rate model for analyzing longitudinal tumor burden (TB) and time-to-event data, as presented in the paper *"Change-Point Detection Using a Cure Rate Joint Model for Longitudinal Tumor Burden and Time-to-Event Data: A Case Study in Modeling Tumor Dynamics in NSCLC Clinical Trial"* by Qu et al. Tailored for oncology clinical trials, such as those in non-small cell lung cancer (NSCLC), this package addresses the dynamic patterns of TB, a critical biomarker reflecting treatment effects. It distinguishes between two patient groups: a "change-point group," where TB initially decreases and later increases (indicating disease progression), and a "stable group," where TB exhibits a sustained decline, potentially representing a cured subset. Leveraging the Monte Carlo Expectation-Maximization (MCEM) algorithm with C++ acceleration, the package efficiently estimates model parameters while incorporating time-to-event data to constrain individual-specific change points in the change-point group, even under censoring. Key functions include `infer_CRJoint_MLE` for parameter estimation and `E_y` for computing expected longitudinal TB trajectories. This flexible and robust framework outperforms traditional models by capturing heterogeneous TB dynamics and providing reliable marginal TB estimates, making it a valuable tool for assessing treatment efficacy in clinical research. The paper link will be available soon.

## Installation

To install the `JoCuR` package from GitHub, use the following commands in R:

```R
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install JoCuR from GitHub (replace with your repository URL)
devtools::install_github("quyixiang/JoCuR")
```


## Usage

### Main Function: `infer_CRJoint_MLE`

This function estimates the parameters of the joint cure rate model using the EM algorithm, accelerated by C++.

#### Example

```R
# Load required library
library(JoCuR)

# Set seed for reproducibility
set.seed(123)

# Simulate example data
sim_data <- CRsimulation(
  fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE),
  fmla.long = as.formula(PCHG ~ 0 + Y0SCALE),
  beta.tte = c(0.2), normal.tte = TRUE, sd.tte = 0.2,
  beta.y = c(-0.5), sd.y = 0.1, beta.y.cure = c(-0.2), sd.y.cured = 0.2,
  cured.rate = 0.4, cured.mean = c(0, -0.2), cured.sd = c(0.2, 0.2),
  cured.corr = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
  randeff.mean = c(0.5, 0, -0.5, 0.5), randeff.sd = rep(0.2, 4),
  randeff.corr = matrix(c(1, -0.4, -0.2, -0.3,
                          -0.4, 1, 0.5, 0.20,
                          -0.2, 0.5, 1, 0.2,
                          -0.3, 0.20, 0.2, 1), nrow = 4),
  n = 100, censor.parameter = 0.5, time.interval = 0.1
)

# Extract datasets
longdat <- sim_data[["longdat"]]
survdat <- sim_data[["survdat"]]

# Define formulas
fmla.tte <- as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE)
fmla.long <- as.formula(PCHG ~ 0 + Y0SCALE)

# Initial values
init <- list(
  mu_r = c(0.5, 0, -0.5, 0.5),
  Sigma_r = diag(rep(0.2, 4)) %*% 
            matrix(c(1, -0.4, -0.2, -0.3,
                     -0.4, 1, 0.5, 0.20,
                     -0.2, 0.5, 1, 0.2,
                     -0.3, 0.20, 0.2, 1), nrow = 4) %*% 
            diag(rep(0.2, 4)),
  beta_y = c(-0.5),
  sigma_y_sq = 0.01,
  beta_tte = c(0.2),
  sigma_tte_sq = 0.04,
  pi_c = 0.4,
  beta_cure = c(-0.2),
  mu_r_cure = c(0, -0.2),
  Sigma_r_cure = diag(c(0.2, 0.2)) %*% 
                 matrix(c(1, 0.5, 0.5, 1), nrow = 2) %*% 
                 diag(c(0.2, 0.2)),
  sigma_y_cure_sq = 0.04
)

# Fit the joint model
results <- infer_CRJoint_MLE(
  survdat = survdat,
  longdat = longdat,
  init = init,
  fmla.tte = fmla.tte,
  fmla.long = fmla.long,
  longdat.time = "visittime",
  id.indicator = "id",
  maxIter = 50,
  tol = 5e-3,
  no_cure = FALSE
)
```

### Additional Function: `E_y`

This function calculates the expected value of the longitudinal outcome at a new time point based on the fitted model.

#### Example

```R
# Assuming 'results' is from the previous example
X_bar <- mean(survdat$Y0SCALE)
W_bar <- mean(survdat$Y0SCALE)
expected_value <- E_y(X_bar, W_bar, results, t_new = 0.1, J = 100)
print(expected_value)
```
