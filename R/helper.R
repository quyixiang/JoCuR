Sigma2P <- function(Sigma) {
  Q <- t(chol(Sigma))
  P <- Q
  diag(P) <- log(diag(Q))
  return(P)
}

vechP2Sigma <- function(vechP) {
  P <- invvech(vechP)
  P[upper.tri(P)] <- 0
  Q <- P2Q(P)
  Sigma <- Q %*% t(Q)
}

Sigma2vechP <- function(Sigma) {
  Q <- t(chol(Sigma))
  P <- Q
  diag(P) <- log(diag(Q))
  return(vech(P))
}

P2Q <- function(P) {
  Q <- P
  diag(Q) <- exp(diag(P))
  return(Q)
}

Z_ci <- function(visittime_vec_i) {
  sec_col <- visittime_vec_i

  result <- as.matrix(cbind(rep(1, length(sec_col)), sec_col))
  colnames(result) <- NULL
  return(as.matrix(result))
}

mu_r_cure_update <- function(ncen, Sigma_r_cure, E_b_cure, E_Delta_cen) {
  first <- sum(E_Delta_cen) * solve(Sigma_r_cure)
  second <- c(0, 0)
  for (i in c(1:ncen)) {
    second <- second + E_Delta_cen[i] * E_b_cure[(2 * i - 1):(2 * i)]
  }
  second <- solve(Sigma_r_cure) %*% second
  return(solve(first) %*% second)
}

Sigma_r_cure_update <- function(ncen, E_b_b_T_mu_cure, E_Delta_cen) {
  first <- matrix(0, 2, 2)
  second <- sum(E_Delta_cen)
  for (i in c(1:ncen)) {
    first <- first + E_Delta_cen[i] * E_b_b_T_mu_cure[(2 * i - 1):(2 * i), ]
  }
  return(first / second)
}

beta_cure_update <- function(ncen, Xcen, ycen, visittime_cen, new_id_cen, E_b_cure, E_Delta_cen) {
  p <- ncol(Xcen)
  a <- matrix(0, p, p)
  b <- matrix(0, p, 1)
  for (i in c(1:ncen)) {
    tmp_visittime <- visittime_cen[new_id_cen == i]
    Z_i_cure <- Z_ci(tmp_visittime)
    x_i <- Xcen[new_id_cen == i, , drop = FALSE]
    y_i <- ycen[new_id_cen == i]

    a <- a + t(x_i) %*% x_i * E_Delta_cen[i]
    b <- b + t(x_i) %*% (y_i - Z_i_cure %*% E_b_cure[(2 * i - 1):(2 * i)]) * E_Delta_cen[i]
  }

  return(solve(a) %*% b)
}

sigma_y_cure_sq_update <- function(ncen, Xcen, ycen, visittime_cen, new_id_cen, beta_cure, mu_r_cure, E_b_cure, E_b_b_T_mu_cure, E_Delta_cen) {
  a <- 0
  b <- 0
  for (i in c(1:ncen)) {
    tmp_visittime <- visittime_cen[new_id_cen == i]
    Z_i_cure <- Z_ci(tmp_visittime)
    E_Zi_bi_value <- Z_i_cure %*% E_b_cure[(2 * i - 1):(2 * i)]

    E_b_b_T_i <- E_b_b_T_mu_cure[(2 * i - 1):(2 * i), ] + 2 * E_b_cure[(2 * i - 1):(2 * i)] %*% t(mu_r_cure) - mu_r_cure %*% t(mu_r_cure)
    x_i <- Xcen[new_id_cen == i, , drop = FALSE]
    y_i <- ycen[new_id_cen == i]
    m_i <- y_i - x_i %*% matrix(beta_cure, ncol = 1)

    a <- a + E_Delta_cen[i] * (tr(m_i %*% t(m_i) - m_i %*% t(E_Zi_bi_value) - E_Zi_bi_value %*% t(m_i)) + tr(t(Z_i_cure) %*% Z_i_cure %*% E_b_b_T_i))
    b <- b + E_Delta_cen[i] * length(y_i)
  }
  return(a / b)
}


wrapped_gradient_function_sigma <- function(par,
                                            Nobs, nobs, E_r_obs, E_r_rT_obs, E_r_rT_mu_obs, Estep_2_obs, yobs, Xobs, Xtte_obs, visittime_obs, tobs, new_id_obs,
                                            Ncen, ncen, E_r_cen, E_r_rT_cen, E_r_rT_mu_cen, Estep_2_cen, ycen, Xcen, Xtte_cen, visittime_cen, E_ti, E_ti_sq, E_g0_ti, E_g1_ti, E_g2_ti, new_id_cen,
                                            beta_tte, sigma_tte_sq, mu_r, beta_y, sigma_y_sq,
                                            E_Delta_cen) {
  Sigma_r <- vechP2Sigma(par)

  obs <- dQdvechP_obs_cpp(Sigma_r, E_r_rT_mu_obs, mu_r, tobs, new_id_obs)
  cen <- dQdvechP_cen_cpp(Sigma_r, E_r_rT_mu_cen, mu_r, E_g2_ti, new_id_cen, E_Delta_cen)
  return(obs + cen)
}


wrapped_gradient_function <- function(par,
                                      Nobs, nobs, E_r_obs, E_r_rT_obs, Estep_2_obs, yobs, Xobs, Xtte_obs, visittime_obs, tobs, new_id_obs,
                                      Ncen, ncen, E_r_cen, E_r_rT_cen, Estep_2_cen, ycen, Xcen, Xtte_cen, visittime_cen, E_ti, E_ti_sq, E_g0_ti, E_g1_ti, E_g2_ti, new_id_cen,
                                      beta_tte, sigma_tte_sq, Sigma_r, beta_y, sigma_y_sq,
                                      E_Delta_cen) {
  obs <- dQdmur_obs_cpp(Sigma_r, E_r_obs, par, tobs, new_id_obs)
  cen <- dQdmur_cen_cpp(Sigma_r, E_r_cen, par, E_g1_ti, new_id_cen, E_Delta_cen)
  return(obs + cen)
}

beta_update <- function(Xobs, Xcen, yobs, ycen, new_id_obs, new_id_cen, E_Z_b_obs, E_Z_b_cen, E_Delta_cen) {
  p <- ncol(Xobs)
  a <- matrix(0, p, p)
  b <- matrix(0, p, 1)
  for (i in c(1:max(new_id_obs))) {
    x_i <- Xobs[new_id_obs == i, , drop = FALSE]
    y_i <- yobs[new_id_obs == i]
    E_Zi_bi_value <- E_Z_b_obs[new_id_obs == i]
    a <- a + t(x_i) %*% x_i
    b <- b + t(x_i) %*% (y_i - E_Zi_bi_value)
  }

  for (i in c(1:max(new_id_cen))) {
    E_1_Delta <- 1 - E_Delta_cen[i]
    x_i <- Xcen[new_id_cen == i, , drop = FALSE]
    y_i <- ycen[new_id_cen == i]
    E_Zi_bi_value <- E_Z_b_cen[new_id_cen == i]
    a <- a + t(x_i) %*% x_i * E_1_Delta
    b <- b + t(x_i) %*% (y_i - E_Zi_bi_value) * E_1_Delta
  }

  return(solve(a) %*% b)
}

sigma_y_sq_update <- function(Estep_2_obs, Estep_2_cen, Nobs, id_number_cen, E_Delta_cen) {
  return(abs((sum(Estep_2_obs) + sum(Estep_2_cen * (1 - E_Delta_cen))) / (Nobs + sum(id_number_cen * (1 - E_Delta_cen)))))
}

beta_tte_update <- function(nobs, ncen, Xtte_obs, Xtte_cen, tobs, E_ti, E_Delta_cen) {
  p <- ncol(Xtte_obs)
  a <- matrix(0, p, p)
  b <- matrix(0, p, 1)
  for (i in 1:nobs) {
    a <- a + (Xtte_obs[i, ]) %*% t(Xtte_obs[i, ])
    b <- b + (Xtte_obs[i, ]) * tobs[i]
  }

  for (i in 1:ncen) {
    E_1_Delta <- 1 - E_Delta_cen[i]
    a <- a + (Xtte_cen[i, ]) %*% t(Xtte_cen[i, ]) * E_1_Delta
    b <- b + (Xtte_cen[i, ]) * E_ti[i] * E_1_Delta
  }
  return(solve(a) %*% b)
}

sigma_tte_sq_update <- function(nobs, ncen, Xtte_obs, Xtte_cen, tobs, E_ti, E_ti_sq, beta_tte, E_Delta_cen) {
  res <- 0

  for (i in 1:nobs) {
    tmp_mu <- t(Xtte_obs[i, ]) %*% beta_tte
    res <- res + (tobs[i] - tmp_mu)^2
  }

  for (i in 1:ncen) {
    E_1_Delta <- 1 - E_Delta_cen[i]
    tmp_mu <- t(Xtte_cen[i, ]) %*% beta_tte
    res <- res + E_1_Delta * (E_ti_sq[i] - 2 * E_ti[i] * tmp_mu + tmp_mu^2)
  }
  return(abs(res[1, 1] / (nobs + ncen - sum(E_Delta_cen))))
}

wrapped_Q_function <- function(x,
                               Nobs, nobs, E_r_obs, E_r_rT_obs, Estep_2_obs, yobs, Xobs, Xtte_obs, visittime_obs, tobs, new_id_obs,
                               Ncen, ncen, E_r_cen, E_r_rT_cen, Estep_2_cen, ycen, Xcen, Xtte_cen, visittime_cen, E_ti, E_ti_sq, E_g0_ti, E_g1_ti, E_g2_ti, new_id_cen,
                               beta_tte, sigma_tte_sq, Sigma_r, beta_y, sigma_y_sq,
                               E_Delta_cen) {
  Q_function_cpp(
    Nobs, nobs, E_r_obs, E_r_rT_obs, Estep_2_obs, yobs, Xobs, Xtte_obs, visittime_obs, tobs, new_id_obs,
    Ncen, ncen, E_r_cen, E_r_rT_cen, Estep_2_cen, ycen, Xcen, Xtte_cen, visittime_cen, E_g0_ti, E_ti, E_ti_sq, new_id_cen,
    beta_tte, sigma_tte_sq, x, Sigma_r, beta_y, sigma_y_sq,
    E_Delta_cen
  )
}


wrapped_Q_function_sigma <- function(x,
                                     Nobs, nobs, E_r_obs, E_r_rT_obs, E_r_rT_mu_obs, Estep_2_obs, yobs, Xobs, Xtte_obs, visittime_obs, tobs, new_id_obs,
                                     Ncen, ncen, E_r_cen, E_r_rT_cen, E_r_rT_mu_cen, Estep_2_cen, ycen, Xcen, Xtte_cen, visittime_cen, E_ti, E_ti_sq, E_g0_ti, E_g1_ti, E_g2_ti, new_id_cen,
                                     beta_tte, sigma_tte_sq, mu_r, beta_y, sigma_y_sq,
                                     E_Delta_cen) {
  Sigma_r <- vechP2Sigma(x)

  Q_function_cpp(
    Nobs, nobs, E_r_obs, E_r_rT_obs, Estep_2_obs, yobs, Xobs, Xtte_obs, visittime_obs, tobs, new_id_obs,
    Ncen, ncen, E_r_cen, E_r_rT_cen, Estep_2_cen, ycen, Xcen, Xtte_cen, visittime_cen, E_g0_ti, E_ti, E_ti_sq, new_id_cen,
    beta_tte, sigma_tte_sq, mu_r, Sigma_r, beta_y, sigma_y_sq,
    E_Delta_cen
  )
}
