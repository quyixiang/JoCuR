#' @title infer_CRJoint_MLE
#' @description This function estimates the parameters of the joint cure rate model using the EM algorithm.
#' @param survdat A data frame containing the survival data.
#' @param longdat A data frame containing the longitudinal data.
#' @param init A list containing the initial values for the parameters.
#' @param fmla.tte A formula specifying the survival outcome.
#' @param fmla.long A formula specifying the longitudinal outcome.
#' @param longdat.time A string specifying the time variable in the longitudinal data.
#' @param id.indicator A string specifying the ID variable in the data.
#' @param maxIter An integer specifying the maximum number of iterations for the EM algorithm.
#' @param tol A numeric value specifying the convergence tolerance for the EM algorithm.
#' @param no_cure A logical value indicating whether to include the cure rate model.
#' @param bounded_slope A logical value indicating whether to use bounded slopes for the longitudinal outcome.
#' @param non_std Integer flag for legacy path: use `0` (default) for standardized identifiable implementation, and `1` for original non-standard implementation.
#' @return A list containing the estimated parameters of the joint cure rate model.
#' @examples
#'   # Load required library
#'   library(JoCuR)
#'
#'   # Set seed for reproducibility
#'   set.seed(123)
#'
#'   # Simulate example data
#'   sim_data <- CRsimulation(
#'     fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE),
#'     fmla.long = as.formula(PCHG ~ 0 + Y0SCALE),
#'     beta.tte = c(0.2), normal.tte = TRUE, sd.tte = 0.2,
#'     beta.y = c(-0.5), sd.y = 0.1, beta.y.cure = c(-0.2), sd.y.cured = 0.2,
#'     cured.rate = 0.4, cured.mean = c(0, -0.2), cured.sd = c(0.2, 0.2),
#'     cured.corr = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
#'     randeff.mean = c(0.5, 0, -0.5, 0.5), randeff.sd = rep(0.2, 4),
#'     randeff.corr = matrix(c(1, -0.4, -0.2, -0.3,
#'                            -0.4, 1, 0.5, 0.20,
#'                            -0.2, 0.5, 1, 0.2,
#'                            -0.3, 0.20, 0.2, 1), nrow = 4),
#'     n = 100, censor.parameter = 0.5, time.interval = 0.1
#'   )
#'
#'   # Extract datasets
#'   longdat <- sim_data[["longdat"]]
#'   survdat <- sim_data[["survdat"]]
#'
#'   # Define formulas
#'   fmla.tte <- as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE)
#'   fmla.long <- as.formula(PCHG ~ 0 + Y0SCALE)
#'
#'   # Initial values
#'   init <- list(
#'     mu_r = c(0.5, 0, -0.5, 0.5),
#'     Sigma_r = diag(rep(0.2, 4)) %*%
#'               matrix(c(1, -0.4, -0.2, -0.3,
#'                       -0.4, 1, 0.5, 0.20,
#'                       -0.2, 0.5, 1, 0.2,
#'                       -0.3, 0.20, 0.2, 1), nrow = 4) %*%
#'               diag(rep(0.2, 4)),
#'     beta_y = c(-0.5),
#'     sigma_y_sq = 0.01,
#'     beta_tte = c(0.2),
#'     sigma_tte_sq = 0.04,
#'     pi_c = 0.4,
#'     beta_cure = c(-0.2),
#'     mu_r_cure = c(0, -0.2),
#'     Sigma_r_cure = diag(c(0.2, 0.2)) %*%
#'                    matrix(c(1, 0.5, 0.5, 1), nrow = 2) %*%
#'                    diag(c(0.2, 0.2)),
#'     sigma_y_cure_sq = 0.04
#'   )
#'
#'   # Fit the joint model
#'   results <- infer_CRJoint_MLE(
#'     survdat = survdat,
#'     longdat = longdat,
#'     init = init,
#'     fmla.tte = fmla.tte,
#'     fmla.long = fmla.long,
#'     longdat.time = "visittime",
#'     id.indicator = "id",
#'     maxIter = 50,
#'     tol = 5e-3,
#'     no_cure = FALSE
#'   )
infer_CRJoint_MLE <- function(survdat, longdat, init, fmla.tte, fmla.long, longdat.time, id.indicator, maxIter = 5000, tol = 8e-4, no_cure = FALSE, bounded_slope = FALSE, non_std = 0) {
  data.list <- load_onearm_data(
    survdat = survdat,
    longdat = longdat,
    fmla.tte = fmla.tte,
    fmla.long = fmla.long,
    longdat.time = longdat.time,
    id.indicator = id.indicator
  )

  data.list[["tcen"]] <- log(data.list[["tcen"]])
  data.list[["tobs"]] <- log(data.list[["tobs"]])

  if (non_std == 1) {
    cure_res <- MCEM_cureJoint(
      data.list,
      initial = init,
      maxIter = maxIter,
      tol = tol,
      no_cure = no_cure,
      bounded_slope = bounded_slope
    )
  } else {
    Xobs <- data.list[["Xlong_obs"]]
    visitobs <- data.list[["visitobs"]]

    Xcen <- data.list[["Xlong_cen"]]
    visitcen <- data.list[["visitcen"]]

    data.list[["Xlong_aug_obs"]] <- cbind(Xobs, Intercept = 1, Time = visitobs)
    data.list[["Xlong_aug_cen"]] <- cbind(Xcen, Intercept = 1, Time = visitcen)

    cure_res <- MCEM_cureJoint_std(
      data.list,
      initial = init,
      maxIter = maxIter,
      tol = tol,
      no_cure = no_cure,
      bounded_slope = bounded_slope
    )
  }

  return(cure_res)
}

infer_CRJoint_MLE_std <- function(survdat, longdat, init, fmla.tte, fmla.long, longdat.time, id.indicator, maxIter = 5000, tol = 8e-4, no_cure = FALSE, bounded_slope = FALSE) {
  infer_CRJoint_MLE(
    survdat = survdat,
    longdat = longdat,
    init = init,
    fmla.tte = fmla.tte,
    fmla.long = fmla.long,
    longdat.time = longdat.time,
    id.indicator = id.indicator,
    maxIter = maxIter,
    tol = tol,
    no_cure = no_cure,
    bounded_slope = bounded_slope,
    non_std = 0
  )
}

#' @title E_y
#' @description This function calculates the expected value of the longitudinal outcome given the cure rate model parameters for a new time point.
#' @param X A vector of covariates for the longitudinal outcome.
#' @param W A vector of covariates for the survival outcome.
#' @param cure_res A list containing the estimated parameters of the cure rate model.
#' @param t_new A new time point for which the expected value is calculated.
#' @param J The number of Monte Carlo samples used for the calculation.
#' @return The expected value of the longitudinal outcome given the cure rate model parameters.
#' @examples
#' \dontrun{
#'  # Make sure that you have run the examples from infer_CRJoint_MLE
#'  X_bar <- mean(survdat$Y0SCALE)
#'  W_bar <- mean(survdat$Y0SCALE)
#'  E_y(X_bar, W_bar, results, 0.1, 100)
#' }
E_y <- function(X, W, cure_res, t_new, J){
  X_bar <- mean(X)
  W_bar <- mean(W)

  idx_w <- 1
  idx_b0 <- 2
  idx_b1 <- 3
  idx_b2 <- 4

  mu_all <- cure_res$mu_r
  Sigma <- cure_res$Sigma_r
  mu_w <- mu_all[idx_w]
  var_w <- Sigma[idx_w, idx_w]

  coef_0 <- Sigma[idx_b0, idx_w] / var_w
  coef_1 <- Sigma[idx_b1, idx_w] / var_w
  coef_2 <- Sigma[idx_b2, idx_w] / var_w

  const_0 <- mu_all[idx_b0] - coef_0 * mu_w
  const_1 <- mu_all[idx_b1] - coef_1 * mu_w
  const_2 <- mu_all[idx_b2] - coef_2 * mu_w

  sum_E0 <- 0
  sum_E1 <- 0
  sum_E2 <- 0

  for (j in 1:J) {
    upper_log <- rnorm(1, mean = cure_res$beta_tte * W_bar, sd = sqrt(cure_res$sigma_tte_sq))
    upper <- exp(upper_log)

    w_samples <- TruncatedNormal::rtnorm(n = 100, mu = mu_w, sd = sqrt(var_w), lb = 0, ub = upper)

    b0_hat <- const_0 + coef_0 * w_samples
    b1_hat <- const_1 + coef_1 * w_samples
    b2_hat <- const_2 + coef_2 * w_samples

    sum_E0 <- sum_E0 + mean(b0_hat)
    sum_E1 <- sum_E1 + mean(b1_hat * (t_new - w_samples) * (t_new <= w_samples))
    sum_E2 <- sum_E2 + mean(b2_hat * (t_new - w_samples) * (t_new > w_samples))
  }

  mean_E0 <- sum_E0 / J
  mean_E1 <- sum_E1 / J
  mean_E2 <- sum_E2 / J

  cp_part <- (X_bar * cure_res$beta_y + mean_E0 + mean_E1 + mean_E2) * (1 - cure_res$pi_c)

  if (length(cure_res$beta_cure) >= 3) {
    stable_part <- (X_bar * cure_res$beta_cure[1] + cure_res$beta_cure[2] + cure_res$beta_cure[3] * t_new) * cure_res$pi_c
  } else {
    stable_part <- (X_bar * cure_res$beta_cure + cure_res$mu_r_cure[1] + cure_res$mu_r_cure[2] * t_new) * cure_res$pi_c
  }

  return(cp_part + stable_part)
}