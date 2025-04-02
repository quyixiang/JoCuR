# library(TruncatedNormal)
# library(randcorr)
# library(dplyr)
# library(Rlab)

rptmvn <- function(mu, sigma, a, b) {
  omegaMean <- mu[1]
  omegaVar <- sigma[1, 1]
  omegaSD <- sqrt(omegaVar)
  ## Check dimensions
  na <- length(a)
  nb <- length(b)
  if (na == 1 & nb > 1) {
    a <- rep(a, nb)
  } else if (na > 1 & nb == 1) {
    b <- rep(b, na)
  } else if (na != nb) {
    stop("dimensions of a and b are not compatible")
  }
  n <- max(na, nb)
  ## Draw omega based on its marginal distribution
  omega <- rtnorm(1, omegaMean, omegaSD, a, b)
  ## Draw b based on its conditional distribution
  bCondCov <- sigma[-1, -1] - tcrossprod(sigma[-1, 1]) / omegaVar
  bCondMean <- t(mu[-1] + tcrossprod(sigma[-1, 1] / sigma[1, 1], omega - omegaMean))
  b <- bCondMean + mvtnorm::rmvnorm(n, sigma = bCondCov)
  colnames(b) <- paste0("b", 1:ncol(b))
  cbind("omega" = omega, "b" = b)
}

rweibullph <- function(X, beta, shape, scale) {
  eta <- (X %*% beta)[, 1]
  logscale <- eta + log(scale)
  scale <- exp(logscale)
  rexp(nrow(X), rate = scale)^(1 / shape)
}


#' Generate Simulated Data for Time-to-Event and Longitudinal Analyses
#'
#' The generate_simulation_data function provides a convenient way to produce simulated data
#' tailored for both time-to-event and longitudinal data analyses. This can be particularly useful for methodological studies,
#' preliminary analyses, and educational purposes.
#'
#' @section Parameters:
#' @param fmla.tte A formula object that delineates the structure of the time-to-event model.
#' @param fmla.long Another formula object that outlines the structure of the longitudinal data model.
#' @param bootstrap A logical value indicating whether to bootstrap the data. If TRUE, the function will sample from the dataset specified in bootstrapfrom. Default is FALSE.
#' @param bootstrapfrom A data frame. When bootstrap is TRUE, this data frame is used for resampling. Default is NULL.
#' @param id.name A character string specifying the name of the ID column in the bootstrapfrom dataset. Default is "id".
#' @param beta.tte A numeric vector that defines the beta coefficients for the TTE model.
#' @param scale.tte A single numeric value that sets the scale parameter for the TTE model.
#' @param shape.tte A single numeric value that sets the shape parameter for the TTE model.
#' @param beta.y A numeric vector that defines the beta coefficients for the longitudinal data model.
#' @param sd.y A single numeric value, representing the standard deviation of the Y variable in the longitudinal data model.
#' @param randeff.mean A numeric vector that defines the average of the random effects.
#' @param randeff.sd A numeric vector that designates the standard deviations associated with the random effects.
#' @param randeff.corr A matrix that sets the correlation structure of the random effects. If NULL, a random structure will be generated.
#' @param n An integer that represents the number of observations to simulate.
#' @param censor.parameter A numeric value used to designate the parameter for the exponential distribution, which is employed to derive censoring times.
#' @param normal.tte Logical indicating if the time-to-event model should follow a log normal distribution. Default is FALSE.
#' @param sd.tte Optional numeric value representing the standard deviation for the time-to-event model. Only used if normal.tte is TRUE.
#' @param time.interval A numeric value that stipulates the expected time interval for visits in the longitudinal data.
#' @param time.interval.sd A numeric value representing the standard deviation for the time interval of visits in the longitudinal data. Default is 0.02.
#' @param seed An integer used to set the seed for random number generation, ensuring reproducibility.
#'
#' @section Returns:
#' The function returns a list comprising three elements:
#' \itemize{
#' \item \code{survdat}: A data frame containing the simulated survival data.
#' \item \code{longdat}: A data frame encompassing the simulated longitudinal data.
#' \item \code{simulation.para}: A list that consolidates the parameters applied during the simulation.
#' }
#'
#' @section Examples:
#' \preformatted{
#' simulation.list <- generate_simulation_data(
#' fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ Y0SCALE + Y1SCALE),
#' fmla.long = as.formula(PCHG ~ 0 + Y0SCALE + Y2SCALE + Y3),
#'
#' beta.tte = c(0.1, 0.05, 0.1), scale.tte = 2, shape.tte = 2,
#' beta.y = c(0.02, -0.02, 0.03), sd.y = 0.1,
#' randeff.mean = c(0.5, 0, -1, 1), randeff.sd = rep(0.2, 4), randeff.corr = NULL,
#' n = 1000, censor.parameter = 2, normal.tte = FALSE, time.interval = 0.1
#' )
#' }
#' @export
CRsimulation <- function(
    fmla.tte, fmla.long,
    bootstrap = FALSE, bootstrapfrom = NULL, id.name = "id",
    beta.tte, scale.tte = NULL, shape.tte = NULL,
    normal.tte = FALSE, sd.tte = NULL,
    beta.y, sd.y, beta.y.cure, sd.y.cured,
    randeff.mean = c(0.5, 0, -1, 1), randeff.sd = rep(0.2, 4), randeff.corr = NULL,
    cured.rate = 0.2, cured.mean = c(0, -0.5), cured.sd = c(0.2, 0.2), cured.corr = NULL,
    n = 100, censor.parameter, time.interval = 0.1, time.interval.sd = 0.02,
    seed = 1) {
  set.seed(seed)

  Xtte.name <- all.vars(fmla.tte)[-c(1:2)]
  Xlong.name <- all.vars(fmla.long)[-1]
  Xtte.response.name <- all.vars(fmla.tte)[1]
  Xtte.indicator.name <- all.vars(fmla.tte)[2]
  Xlong.response.name <- all.vars(fmla.long)[1]
  Xall.name <- unique(c(Xtte.name, Xlong.name))
  if (bootstrap) {
    id.sel <- sample(bootstrapfrom[[id.name]], size = n, replace = TRUE)

    X <- data.frame(matrix(ncol = length(Xall.name), nrow = 0))
    colnames(X) <- Xall.name

    for (id in id.sel) {
      X <- rbind(X, bootstrapfrom[bootstrapfrom[[id.name]] == id, Xall.name, drop = FALSE]) %>% as.matrix()
    }
    row.names(X) <- c(1:n)
  } else {
    X <- matrix(rnorm(n * length(Xall.name), mean = 0, sd = 1), nrow = n, ncol = length(Xall.name), byrow = FALSE)
  }
  colnames(X) <- Xall.name
  X <- data.frame(X)
  fmla.tte.rhs <- as.formula(paste("~", deparse(fmla.tte[[3]])))
  fmla.long.rhs <- as.formula(paste("~", deparse(fmla.long[[3]])))
  Xtte <- model.matrix(fmla.tte.rhs, data = X)
  Xlong <- model.matrix(fmla.long.rhs, data = X)
  Xlong.name.intercept <- colnames(Xlong)
  X_combined <- cbind(Xtte, Xlong)
  X_combined_df <- as.data.frame(X_combined)
  X_combined_unique <- X_combined_df %>%
    select(unique(colnames(X_combined_df)))
  X.model <- as.matrix(X_combined_unique)

  id <- c(1:n)
  survdat <- as.data.frame(X.model)
  survdat$id <- id
  survdat <- survdat[, c("id", names(survdat)[!names(survdat) %in% "id"])]
  if (normal.tte) {
    survdat$event_years <- exp(rnorm(n = n, mean = Xtte %*% beta.tte, sd = sd.tte))
  } else {
    survdat$event_years <- rweibullph(Xtte, beta.tte, scale.tte, scale.tte)
  }

  survdat$id <- 1:nrow(survdat)

  survdat$cured <- rbern(n, cured.rate)
  survdat[survdat$cured == 1, "event_years"] <- Inf

  set.seed(seed)
  survdat$censor_years <- rexp(nrow(survdat), censor.parameter)
  survdat[[Xtte.indicator.name]] <- ifelse(survdat$event_years <= survdat$censor_years, 1, 0)
  survdat[[Xtte.response.name]] <- pmin(survdat$censor_years, survdat$event_years)


  if (!is.null(randeff.corr)) {
    randeff.corr <- matrix(randeff.corr, nrow = 4, ncol = 4)
  } else {
    set.seed(seed)
    randeff.corr <- randcorr(4)
  }

  randeff.cov <- diag(randeff.sd) %*% randeff.corr %*% diag(randeff.sd)
  randeff <- rptmvn(randeff.mean, randeff.cov, 0, survdat$event_years)

  if (!is.null(cured.corr)) {
    cured.corr <- matrix(cured.corr, nrow = 2, ncol = 2)
  } else {
    set.seed(seed)
    cured.corr <- randcorr(2)
  }

  cured.cov <- diag(cured.sd) %*% cured.corr %*% diag(cured.sd)
  randeff.cured <- mvtnorm::rmvnorm(n, mean = cured.mean, sigma = cured.cov)
  colnames(randeff.cured) <- c("b0_cured", "b1_cured")
  survdat <- cbind(survdat, randeff, randeff.cured)

  survdat$nvisits <- ceiling(survdat[[Xtte.response.name]] / time.interval)

  survdat <- survdat %>% filter(nvisits > 0)


  longdat <- survdat %>%
    filter(nvisits > 0)
  rep.indx <- lapply(1:nrow(longdat), function(i) rep(longdat$id[i], longdat$nvisits[i])) %>% unlist()

  set.seed(seed)
  longdat <- longdat[rep.indx, ] %>%
    group_by(id) %>%
    mutate(visitnum = row_number()) %>%
    rowwise() %>%
    mutate(rand_val = abs(rnorm(1, mean = 0, sd = time.interval.sd))) %>%
    mutate(visittime = abs(time.interval * visitnum - rand_val))

  longdat <- longdat %>%
    mutate(visittime = if_else(nvisits %in% c(1, 2) & visitnum == 1,
      0.1 * get(as.character(Xtte.response.name)),
      visittime
    ))

  longdat <- longdat %>%
    ungroup() %>%
    mutate(
      delta = visittime - omega,
      eta_re = b1 + if_else(delta < 0, b2 * delta, b3 * delta)
    ) %>%
    rowwise() %>%
    mutate(
      eta_fe = sum(c_across(all_of(Xlong.name.intercept)) * beta.y),
      y_uncured = rnorm(1, mean = eta_fe + eta_re, sd = sd.y)
    )

  longdat <- longdat %>%
    rowwise() %>%
    mutate(
      eta_fe_cure = sum(c_across(all_of(Xlong.name.intercept)) * beta.y.cure),
      y_cured = eta_fe_cure + rnorm(1, mean = b0_cured + b1_cured * visittime, sd = sd.y.cured)
    )

  longdat <- longdat %>%
    mutate(
      !!sym(Xlong.response.name) := case_when(
        cured == 0 ~ y_uncured,
        TRUE ~ y_cured
      )
    )
  longdat <- longdat %>% filter(visittime <= !!sym(Xtte.response.name))

  longdat <- longdat %>%
    group_by(id) %>%
    mutate(nvisits_aftercensor = max(row_number())) %>%
    ungroup()
  survdat <- merge(survdat, longdat %>% select(id, nvisits_aftercensor) %>% distinct(id, .keep_all = TRUE))


  cat(
    paste0(
      "Number of Censor: ", sum(survdat[[Xtte.indicator.name]] == 0), "\nNumber of Observation: ", sum(survdat[[Xtte.indicator.name]] == 1),
      "\nProportion of Censor: ", sum(survdat[[Xtte.indicator.name]] == 0) / nrow(survdat), "\n"
    )
  )

  rownames(randeff.corr) <- c("omega", "b1", "b2", "b3")
  colnames(randeff.corr) <- c("omega", "b1", "b2", "b3")


  simulation.para <- c(
    beta.tte = beta.tte,
    scale.tte = scale.tte,
    shape.tte = shape.tte,
    beta.y = beta.y,
    sd.y = sd.y,
    randeff.mean = randeff.mean,
    randeff.sd = randeff.sd,
    randeff.corr = randeff.corr,
    cured.rate = cured.rate
  )
  return(list(survdat = survdat, longdat = longdat, simulation.para = simulation.para))
}



CRsimulation_fixed_visittime <- function(
    fmla.tte, fmla.long,
    bootstrap = FALSE, bootstrapfrom = NULL, id.name = "id",
    beta.tte, scale.tte = NULL, shape.tte = NULL,
    normal.tte = FALSE, sd.tte = NULL,
    beta.y, sd.y, beta.y.cure, sd.y.cured,
    randeff.mean = c(0.5, 0, -1, 1), randeff.sd = rep(0.2, 4), randeff.corr = NULL,
    cured.rate = 0.2, cured.mean = c(0, -0.5), cured.sd = c(0.2, 0.2), cured.corr = NULL,
    n = 100, censor.parameter, time.interval = 0.1, time.interval.sd = 0.02, max_visittime = 2,
    seed = 1) {
  set.seed(seed)

  Xtte.name <- all.vars(fmla.tte)[-c(1:2)]
  Xlong.name <- all.vars(fmla.long)[-1]
  Xtte.response.name <- all.vars(fmla.tte)[1]
  Xtte.indicator.name <- all.vars(fmla.tte)[2]
  Xlong.response.name <- all.vars(fmla.long)[1]
  Xall.name <- unique(c(Xtte.name, Xlong.name))
  if (bootstrap) {
    id.sel <- sample(bootstrapfrom[[id.name]], size = n, replace = TRUE)

    X <- data.frame(matrix(ncol = length(Xall.name), nrow = 0))
    colnames(X) <- Xall.name

    for (id in id.sel) {
      X <- rbind(X, bootstrapfrom[bootstrapfrom[[id.name]] == id, Xall.name, drop = FALSE]) %>% as.matrix()
    }
    row.names(X) <- c(1:n)
  } else {
    X <- matrix(rnorm(n * length(Xall.name), mean = 0, sd = 1), nrow = n, ncol = length(Xall.name), byrow = FALSE)
  }
  colnames(X) <- Xall.name
  X <- data.frame(X)
  fmla.tte.rhs <- as.formula(paste("~", deparse(fmla.tte[[3]])))
  fmla.long.rhs <- as.formula(paste("~", deparse(fmla.long[[3]])))
  Xtte <- model.matrix(fmla.tte.rhs, data = X)
  Xlong <- model.matrix(fmla.long.rhs, data = X)
  Xlong.name.intercept <- colnames(Xlong)
  X_combined <- cbind(Xtte, Xlong)
  X_combined_df <- as.data.frame(X_combined)
  X_combined_unique <- X_combined_df %>%
    select(unique(colnames(X_combined_df)))
  X.model <- as.matrix(X_combined_unique)

  id <- c(1:n)
  survdat <- as.data.frame(X.model)
  survdat$id <- id
  survdat <- survdat[, c("id", names(survdat)[!names(survdat) %in% "id"])]
  if (normal.tte) {
    survdat$event_years <- exp(rnorm(n = n, mean = Xtte %*% beta.tte, sd = sd.tte))
  } else {
    survdat$event_years <- rweibullph(Xtte, beta.tte, scale.tte, scale.tte)
  }

  survdat$id <- 1:nrow(survdat)

  survdat$cured <- rbern(n, cured.rate)
  survdat[survdat$cured == 1, 'event_years'] <- Inf

  set.seed(seed)
  survdat$censor_years <- rexp(nrow(survdat), censor.parameter)
  survdat[[Xtte.indicator.name]] <- ifelse(survdat$event_years <= survdat$censor_years, 1, 0)
  survdat[[Xtte.response.name]] <- pmin(survdat$censor_years, survdat$event_years)


  if (!is.null(randeff.corr)) {
    randeff.corr <- matrix(randeff.corr, nrow = 4, ncol = 4)
  } else {
    set.seed(seed)
    randeff.corr <- randcorr(4)
  }

  randeff.cov <- diag(randeff.sd) %*% randeff.corr %*% diag(randeff.sd)
  randeff <- rptmvn(randeff.mean, randeff.cov, 0, survdat$event_years)

  if (!is.null(cured.corr)){
    cured.corr <- matrix(cured.corr, nrow = 2, ncol = 2)
  } else {
    set.seed(seed)
    cured.corr <- randcorr(2)
  }

  cured.cov <- diag(cured.sd) %*% cured.corr %*% diag(cured.sd)
  randeff.cured <- mvtnorm::rmvnorm(n, mean = cured.mean, sigma = cured.cov)
  colnames(randeff.cured) <- c("b0_cured", "b1_cured")
  survdat <- cbind(survdat, randeff, randeff.cured)

  survdat$nvisits <- ceiling(max_visittime / time.interval)

  survdat <- survdat %>% filter(nvisits > 0)


  longdat <- survdat %>%
    filter(nvisits > 0)
  rep.indx <- lapply(1:nrow(longdat), function(i) rep(longdat$id[i], longdat$nvisits[i])) %>% unlist()

  set.seed(seed)
  longdat <- longdat[rep.indx, ] %>%
    group_by(id) %>%
    mutate(visitnum = row_number()) %>%
    rowwise() %>%
    mutate(rand_val = abs(rnorm(1, mean = 0, sd = time.interval.sd))) %>%
    mutate(visittime = abs(time.interval * visitnum))

  longdat <- longdat %>%
    mutate(visittime = if_else(nvisits %in% c(1,2) & visitnum == 1,
                               0.1 * get(as.character(Xtte.response.name)),
                               visittime))

  longdat <- longdat %>%
    ungroup() %>%
    mutate(
      delta = visittime - omega,
      eta_re = b1 + if_else(delta < 0, b2 * delta, b3 * delta)
    ) %>%
    rowwise() %>%
    mutate(
      eta_fe = sum(c_across(all_of(Xlong.name.intercept)) * beta.y),
      y_uncured = rnorm(1, mean = eta_fe + eta_re, sd = sd.y)
    )

  longdat <- longdat %>%
    rowwise() %>%
    mutate(
      eta_fe_cure = sum(c_across(all_of(Xlong.name.intercept)) * beta.y.cure),
      y_cured = eta_fe_cure + rnorm(1, mean = b0_cured + b1_cured * visittime, sd = sd.y.cured))

  longdat <- longdat %>%
    mutate(
      !!sym(Xlong.response.name) := case_when(
        cured == 0 ~ y_uncured,
        TRUE ~ y_cured
      )
    )

  longdat <- longdat %>%
    group_by(id) %>%
    mutate(nvisits_aftercensor = max(row_number())) %>%
    ungroup()
  survdat <- merge(survdat, longdat %>% select(id, nvisits_aftercensor) %>% distinct(id, .keep_all = TRUE))


  cat(
    paste0(
      "Number of Censor: ", sum(survdat[[Xtte.indicator.name]] == 0), "\nNumber of Observation: ", sum(survdat[[Xtte.indicator.name]] == 1),
      "\nProportion of Censor: ", sum(survdat[[Xtte.indicator.name]] == 0) / nrow(survdat), "\n"
    )
  )

  rownames(randeff.corr) <- c("omega", "b1", "b2", "b3")
  colnames(randeff.corr) <- c("omega", "b1", "b2", "b3")


  simulation.para <- c(
    beta.tte = beta.tte,
    scale.tte = scale.tte,
    shape.tte = shape.tte,
    beta.y = beta.y,
    sd.y = sd.y,
    randeff.mean = randeff.mean,
    randeff.sd = randeff.sd,
    randeff.corr = randeff.corr,
    cured.rate = cured.rate
  )
  return(list(survdat = survdat, longdat = longdat, simulation.para = simulation.para))
}
