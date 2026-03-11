testthat::test_that("simulation, estimation, and E_y smoke test", {
  testthat::skip_on_cran()

  set.seed(123)

  fmla.tte <- as.formula(survival::Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE)
  fmla.long <- as.formula(PCHG ~ 0 + Y0SCALE)

  sim_data <- JoCuR::CRsimulation(
    fmla.tte = fmla.tte,
    fmla.long = fmla.long,
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

  testthat::expect_true(all(c("survdat", "longdat") %in% names(sim_data)))

  survdat <- sim_data[["survdat"]]
  longdat <- sim_data[["longdat"]]

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

  suppressWarnings(capture.output({
    fit <- JoCuR::infer_CRJoint_MLE(
      survdat = survdat,
      longdat = longdat,
      init = init,
      fmla.tte = fmla.tte,
      fmla.long = fmla.long,
      longdat.time = "visittime",
      id.indicator = "id",
      maxIter = 3,
      tol = 1e-2,
      no_cure = FALSE
    )
  }))

  testthat::expect_true(is.list(fit))
  testthat::expect_true(all(c("mu_r", "Sigma_r", "beta_y", "pi_c", "beta_cure") %in% names(fit)))

  ey <- JoCuR::E_y(
    X = survdat$Y0SCALE,
    W = survdat$Y0SCALE,
    cure_res = fit,
    t_new = 0.5,
    J = 30
  )

  testthat::expect_true(is.numeric(ey))
  testthat::expect_length(ey, 1)
  testthat::expect_true(is.finite(ey))
})
