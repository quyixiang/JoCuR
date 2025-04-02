

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
infer_CRJoint_MLE <- function(survdat, longdat, init, fmla.tte, fmla.long, longdat.time, id.indicator, maxIter = 5000, tol = 8e-4, no_cure = FALSE) {
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


  cure_res <- MCEM_cureJoint(data.list, initial = init, maxIter = maxIter, tol = tol, no_cure = no_cure)
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
#'  # Make sure that you have run the examples from infer_CRJoint_MLE
#'  X_bar <- mean(survdat$Y0SCALE)
#'  W_bar <- mean(survdat$Y0SCALE)
#'  E_y(X_bar, W_bar, results, 0.1, 100)

E_y <- function(X, W, cure_res, t_new, J){
  X_bar = mean(X)
  W_bar = mean(W)
  E0 = rep(0, J)
  E1 = rep(0, J)
  E2 = rep(0, J)
  
  for (j in c(1:J)){
    upper_log = rnorm(1, mean = cure_res$beta_tte * W_bar, sd = sqrt(cure_res$sigma_tte_sq))
    upper = exp(upper_log)
    rf = TruncatedNormal::rtmvnorm(100, mu = cure_res$mu_r, sigma = cure_res$Sigma_r, lb = c(0, -Inf, -Inf, -Inf), ub = c(upper, Inf, Inf, Inf))
    E0[j] = mean(rf[,2])
    E1[j] = mean(rf[,3] * (t_new - rf[,1]) * ( t_new <= rf[,1]))
    E2[j] = mean(rf[,4] * (t_new - rf[,1]) * ( t_new > rf[,1]))
    
  }
  res = (X_bar * cure_res$beta_y + mean(E0) + mean(E1) + mean(E2)) * (1-cure_res$pi_c) + (X_bar * cure_res$beta_cure + cure_res$mu_r_cure[1] + cure_res$mu_r_cure[2] * t_new) * cure_res$pi_c
  return(unname(res))
}