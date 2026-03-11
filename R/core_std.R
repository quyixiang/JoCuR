MCEM_cureJoint_std <- function(data.list, tol = 1e-6, maxIter = 1000, initial = NULL, gamma = 0.0001, no_cure = FALSE, bounded_slope = FALSE) {
  nobs <- data.list[["nobs"]]
  ncen <- data.list[["ncen"]]
  visittime_obs <- data.list[["visitobs"]]
  visittime_cen <- data.list[["visitcen"]]
  new_id_obs <- data.list[["idobs"]]
  new_id_cen <- data.list[["idcen"]]
  id_number_cen <- as.vector(table(new_id_cen))
  

  Xobs <- data.list[["Xlong_obs"]]
  Xcen <- data.list[["Xlong_cen"]]
  

  Xobs_aug <- data.list[["Xlong_aug_obs"]]
  Xcen_aug <- data.list[["Xlong_aug_cen"]]
  
  Xtte_obs <- data.list[["Xtte_obs"]]
  Xtte_cen <- data.list[["Xtte_cen"]]
  yobs <- data.list[["yobs"]]
  ycen <- data.list[["ycen"]]
  # this t is in log scale
  tobs <- data.list[["tobs"]]
  tcen <- data.list[["tcen"]]
  Nobs <- data.list[["Nobs"]]
  Ncen <- data.list[["Ncen"]]
  
  p <- dim(Xobs)[2]       # Original dimension
  p_aug <- dim(Xobs_aug)[2] # Augmented dimension (p + 2)
  ptte <- dim(Xtte_obs)[2]

  # --- INITIALIZATION BLOCK ---
  if (is.null(initial)) {
    Sigma_r_Q <- diag(c(0.3, 0.3, 0.3, 0.3))
    Sigma_r <- Sigma_r_Q %*% t(Sigma_r_Q)
    mu_r <- c(0.5, 0, -0.5, 0.5)
    
    beta_y <- rep(0, p)
    sigma_y_sq <- 1
    beta_tte <- rep(0, ptte)
    sigma_tte_sq <- 0.06
    
    # [CHANGE] Random Mean fixed to 0 for identifiability
    mu_r_cure <- rep(0, 2) 
    Sigma_r_cure <- diag(c(0.2, 0.2))
    
    # [CHANGE] beta_cure now has dimension p_aug (includes int + slope)
    beta_cure_aug <- rep(0, p_aug)
    
    sigma_y_cure_sq <- 0.01
    if (!(no_cure)) {
      pi_c <- 0.4
    } else {
      pi_c <- 0
    }
  } else {
    mu_r <- initial[["mu_r"]]
    Sigma_r <- initial[["Sigma_r"]]
    beta_y <- initial[["beta_y"]]
    sigma_y_sq <- initial[["sigma_y_sq"]]
    beta_tte <- initial[["beta_tte"]]
    sigma_tte_sq <- initial[["sigma_tte_sq"]]

    if (!(no_cure)) {
      pi_c <- initial[["pi_c"]]
      
      # [CHANGE] Handle initial beta_cure dimension mismatch
      # If user provided p-dim beta, pad it. If p_aug-dim, keep it.
      init_beta_cure <- initial[["beta_cure"]]
      if(length(init_beta_cure) == p) {
         beta_cure_aug <- c(init_beta_cure, 0, 0)
      } else {
         beta_cure_aug <- init_beta_cure
      }

      # [CHANGE] Force mu_r_cure to 0
      mu_r_cure <- rep(0, 2)
      
      Sigma_r_cure <- initial[["Sigma_r_cure"]]
      sigma_y_cure_sq <- initial[["sigma_y_cure_sq"]]
    } else {
      pi_c <- 0
      mu_r_cure <- rep(0, 2)
      Sigma_r_cure <- diag(c(0.2, 0.2))
      beta_cure_aug <- rep(0, p_aug)
      sigma_y_cure_sq <- 0.01
    }
  }
  
  # Construct beta_y_aug (padded with 0s for the new columns)
  # This ensures beta_y ignores the Intercept/Time columns in C++
  beta_y_aug <- c(beta_y, 0, 0)

  vechP <- vech(Sigma2P(Sigma_r))

  print("Conducting EM algorithm, which may finish early.")
  pb <- utils::txtProgressBar(min = 0, max = maxIter, style = 3)
  iter <- 1
  eps <- Inf
  while (iter <= maxIter && eps > tol) {
    
    # [CHANGE] Pass AUGMENTED matrices (from data.list) and PADDED betas to C++
    # This allows Stable group to see Int/Time, while ChangePoint group ignores them (via 0 coefs)
    
    E_loop_obs_cpp_res <- E_loop_obs_cpp(visittime_obs, new_id_obs, Xobs_aug, yobs, Xtte_obs, tobs, as.vector(beta_y_aug), as.vector(beta_tte), sigma_tte_sq, sigma_y_sq, mu_r, Sigma_r, nobs, Nobs, 100)

    E_r_obs <- E_loop_obs_cpp_res[["E_r_obs"]]
    E_r_rT_obs <- E_loop_obs_cpp_res[["E_r_rT_obs"]]
    E_r_rT_mu_obs <- E_loop_obs_cpp_res[["E_r_rT_mu_obs"]]
    E_Z_b_obs <- E_loop_obs_cpp_res[["E_Z_b_obs"]]
    Estep_2_obs <- E_loop_obs_cpp_res[["Estep_2_obs"]]

    E_loop_cpp_res <- E_loop_cpp(visittime_cen, new_id_cen, Xcen_aug, ycen, Xtte_cen, tcen, as.vector(beta_y_aug), as.vector(beta_tte), sigma_tte_sq, sigma_y_sq, mu_r, Sigma_r, mu_r_cure, Sigma_r_cure, as.vector(beta_cure_aug), sigma_y_cure_sq, ncen, Ncen, 100)

    E_denominator_cen <- E_loop_cpp_res[["E_denominator_cen"]]
    E_Delta1 <- E_loop_cpp_res[["E_Delta1"]]
    E_t_cen <- E_loop_cpp_res[["E_t_cen"]]
    E_t_sq_cen <- E_loop_cpp_res[["E_t_sq_cen"]]
    E_g0_t_cen <- E_loop_cpp_res[["E_g0_t_cen"]]
    E_g1_t_cen <- E_loop_cpp_res[["E_g1_t_cen"]]
    E_g2_t_cen <- E_loop_cpp_res[["E_g2_t_cen"]]
    E_r_cen <- E_loop_cpp_res[["E_r_cen"]]
    E_r_rT_cen <- E_loop_cpp_res[["E_r_rT_cen"]]
    E_r_rT_mu_cen <- E_loop_cpp_res[["E_r_rT_mu_cen"]]
    E_Z_b_cen <- E_loop_cpp_res[["E_Z_b_cen"]]
    Estep_2_cen <- E_loop_cpp_res[["Estep_2_cen"]]
    if (!(no_cure)) {
      E_b_cure <- E_loop_cpp_res[["E_b_cure"]]
      E_b_b_T_mu_cure <- E_loop_cpp_res[["E_b_b_T_mu_cure"]]
    }

    E_Delta_cen <- pi_c * E_Delta1 / (pi_c * E_Delta1 + (1 - pi_c) * E_denominator_cen)

    # M step
    pi_c_new <- sum(E_Delta_cen) / (ncen + nobs)

    # [CHANGE] Update beta_y using ORIGINAL matrices (Xobs, Xcen)
    # This keeps the change-point model structure exactly as it was.
    beta_y_new <- beta_update(Xobs, Xcen, yobs, ycen, new_id_obs, new_id_cen, E_Z_b_obs, E_Z_b_cen, E_Delta_cen)
    
    # [CHANGE] Re-pad beta_y for the next C++ call
    beta_y_aug_new <- c(beta_y_new, 0, 0)
    
    sigma_y_sq_new <- sigma_y_sq_update(Estep_2_obs, Estep_2_cen, Nobs, id_number_cen, E_Delta_cen)

    beta_tte_new <- beta_tte_update(nobs, ncen, Xtte_obs, Xtte_cen, tobs, E_t_cen, E_Delta_cen)
    sigma_tte_sq_new <- sigma_tte_sq_update(nobs, ncen, Xtte_obs, Xtte_cen, tobs, E_t_cen, E_t_sq_cen, beta_tte, E_Delta_cen)

    if (!(no_cure)) {
      mu_r_cure_new <- rep(0, 2) # Fixed at 0
      Sigma_r_cure_new <- Sigma_r_cure
      beta_cure_aug_new <- beta_cure_aug
      sigma_y_cure_sq_new <- sigma_y_cure_sq
      tryCatch(
        {
          # mu_r_cure_new <- mu_r_cure_update(...) # Skipped
          
          Sigma_r_cure_new <- Sigma_r_cure_update(ncen, E_b_b_T_mu_cure, E_Delta_cen)
          
          # [CHANGE] Update beta_cure using AUGMENTED matrix (Xcen_aug)
          # This estimates Covariates + Intercept + Time
          beta_cure_aug_new <- beta_cure_update(ncen, Xcen_aug, ycen, visittime_cen, new_id_cen, E_b_cure, E_Delta_cen)
          
          # [CHANGE] Use Xcen_aug and beta_cure_aug for variance update
          sigma_y_cure_sq_new <- max(sigma_y_cure_sq_update(ncen, Xcen_aug, ycen, visittime_cen, new_id_cen, beta_cure_aug, mu_r_cure, E_b_cure, E_b_b_T_mu_cure, E_Delta_cen), 1e-6)
        },
        error = function(e) {
        }
      )
    } else {
      mu_r_cure_new <- mu_r_cure
      Sigma_r_cure_new <- Sigma_r_cure
      beta_cure_aug_new <- beta_cure_aug
      sigma_y_cure_sq_new <- sigma_y_cure_sq
    }

    # [CHANGE] Q Function uses Augmented matrices and padded betas
    Q_function_value <- wrapped_Q_function(
      x = mu_r,
      Nobs = Nobs, nobs = nobs, E_r_obs = E_r_obs, E_r_rT_obs = E_r_rT_obs, Estep_2_obs = Estep_2_obs, yobs = yobs, Xobs = Xobs_aug, Xtte_obs = Xtte_obs,
      visittime_obs = visittime_obs, tobs = tobs, new_id_obs = new_id_obs,
      Ncen = Ncen, ncen = ncen, E_r_cen = E_r_cen, E_r_rT_cen = E_r_rT_cen, Estep_2_cen = Estep_2_cen, ycen = ycen, Xcen = Xcen_aug, Xtte_cen = Xtte_cen,
      visittime_cen = visittime_cen, E_ti = E_t_cen, E_ti_sq = E_t_sq_cen, E_g0_ti = E_g0_t_cen, E_g1_ti = E_g1_t_cen, E_g2_ti = E_g2_t_cen, new_id_cen = new_id_cen,
      beta_tte = beta_tte, sigma_tte_sq = sigma_tte_sq, Sigma_r = Sigma_r, beta_y = beta_y_aug, sigma_y_sq = sigma_y_sq, E_Delta_cen = E_Delta_cen
    )

    # Optimization blocks passed augmented X and padded betas
    if (!(bounded_slope)) {
      mu_r_new <- mu_r
      tryCatch(
        {
          mu_optim <- optim(
            par = mu_r, fn = wrapped_Q_function, gr = wrapped_gradient_function,
            lower = c(0, -Inf, -Inf, -Inf), method = "L-BFGS-B", control = list(fnscale = -1),
            Nobs = Nobs, nobs = nobs, E_r_obs = E_r_obs, E_r_rT_obs = E_r_rT_obs, Estep_2_obs = Estep_2_obs, yobs = yobs, Xobs = Xobs_aug, Xtte_obs = Xtte_obs,
            visittime_obs = visittime_obs, tobs = tobs, new_id_obs = new_id_obs,
            Ncen = Ncen, ncen = ncen, E_r_cen = E_r_cen, E_r_rT_cen = E_r_rT_cen, Estep_2_cen = Estep_2_cen, ycen = ycen, Xcen = Xcen_aug, Xtte_cen = Xtte_cen,
            visittime_cen = visittime_cen, E_ti = E_t_cen, E_ti_sq = E_t_sq_cen, E_g0_ti = E_g0_t_cen, E_g1_ti = E_g1_t_cen, E_g2_ti = E_g2_t_cen, new_id_cen = new_id_cen,
            beta_tte = beta_tte, sigma_tte_sq = sigma_tte_sq, Sigma_r = Sigma_r, beta_y = beta_y_aug, sigma_y_sq = sigma_y_sq, E_Delta_cen = E_Delta_cen
          )
          mu_r_new <- mu_optim$par
        },
        error = function(e) {
        }
      )
    } else {
       # (Same augmentation for bounded_slope block)
       mu_r_new <- mu_r
       tryCatch(
        {
          mu_optim <- optim(
            par = mu_r, fn = wrapped_Q_function, gr = wrapped_gradient_function,
            lower = c(0, -Inf, -Inf, 0), upper = c(Inf, Inf, 0, Inf),
            method = "L-BFGS-B", control = list(fnscale = -1),
            Nobs = Nobs, nobs = nobs, E_r_obs = E_r_obs, E_r_rT_obs = E_r_rT_obs, Estep_2_obs = Estep_2_obs, yobs = yobs, Xobs = Xobs_aug, Xtte_obs = Xtte_obs,
            visittime_obs = visittime_obs, tobs = tobs, new_id_obs = new_id_obs,
            Ncen = Ncen, ncen = ncen, E_r_cen = E_r_cen, E_r_rT_cen = E_r_rT_cen, Estep_2_cen = Estep_2_cen, ycen = ycen, Xcen = Xcen_aug, Xtte_cen = Xtte_cen,
            visittime_cen = visittime_cen, E_ti = E_t_cen, E_ti_sq = E_t_sq_cen, E_g0_ti = E_g0_t_cen, E_g1_ti = E_g1_t_cen, E_g2_ti = E_g2_t_cen, new_id_cen = new_id_cen,
            beta_tte = beta_tte, sigma_tte_sq = sigma_tte_sq, Sigma_r = Sigma_r, beta_y = beta_y_aug, sigma_y_sq = sigma_y_sq, E_Delta_cen = E_Delta_cen
          )
          mu_r_new <- mu_optim$par
        },
        error = function(e) {
        }
      )
    }

    vechP_new <- vechP
    Sigma_r_new <- vechP2Sigma(vechP)
    tryCatch(
      {
        P_optim <- optim(
          par = vechP, fn = wrapped_Q_function_sigma, gr = wrapped_gradient_function_sigma,
          method = "L-BFGS-B", 
          lower = rep(-2, 10), upper = rep(2, 10),
          control = list(fnscale = -1),
          Nobs = Nobs, nobs = nobs, E_r_obs = E_r_obs, E_r_rT_obs = E_r_rT_obs, E_r_rT_mu_obs = E_r_rT_mu_obs, Estep_2_obs = Estep_2_obs, yobs = yobs, Xobs = Xobs_aug, Xtte_obs = Xtte_obs,
          visittime_obs = visittime_obs, tobs = tobs, new_id_obs = new_id_obs,
          Ncen = Ncen, ncen = ncen, E_r_cen = E_r_cen, E_r_rT_cen = E_r_rT_cen, E_r_rT_mu_cen = E_r_rT_mu_cen, Estep_2_cen = Estep_2_cen, ycen = ycen, Xcen = Xcen_aug, Xtte_cen = Xtte_cen,
          visittime_cen = visittime_cen, E_ti = E_t_cen, E_ti_sq = E_t_sq_cen, E_g0_ti = E_g0_t_cen, E_g1_ti = E_g1_t_cen, E_g2_ti = E_g2_t_cen, new_id_cen = new_id_cen,
          beta_tte = beta_tte, sigma_tte_sq = sigma_tte_sq, beta_y = beta_y_aug, sigma_y_sq = sigma_y_sq, E_Delta_cen = E_Delta_cen,
          mu_r = mu_r
        )

        vechP_new <- P_optim$par

        Sigma_r_new <- vechP2Sigma(vechP_new)
      },
      error = function(e) {
      }
    )

    # Use original beta_y (not augmented) for convergence check and return, to keep format clean
    par_old <- c(mu_r, vech(Sigma_r), beta_y, sigma_y_sq, beta_tte, sigma_tte_sq, pi_c, beta_cure_aug, mu_r_cure, vech(Sigma_r_cure), sigma_y_cure_sq)
    par_new <- c(mu_r_new, vech(Sigma_r_new), beta_y_new, sigma_y_sq_new, beta_tte_new, sigma_tte_sq_new, pi_c_new, beta_cure_aug_new, mu_r_cure_new, vech(Sigma_r_cure_new), sigma_y_cure_sq_new)

    eps <- sqrt(abs(sum((par_old - par_new)^2) / (sum(par_old))^2))
    
    # update
    mu_r <- mu_r_new
    vechP <- vechP_new
    Sigma_r <- Sigma_r_new

    mu_r_cure <- mu_r_cure_new
    Sigma_r_cure <- Sigma_r_cure_new

    beta_y <- beta_y_new
    beta_y_aug <- beta_y_aug_new # Update the padded version for next loop
    
    sigma_y_sq <- sigma_y_sq_new
    
    beta_cure_aug <- beta_cure_aug_new
    sigma_y_cure_sq <- sigma_y_cure_sq_new

    beta_tte <- beta_tte_new
    sigma_tte_sq <- sigma_tte_sq_new
    pi_c <- pi_c_new

    utils::setTxtProgressBar(pb, iter)
    iter <- iter + 1
  }
  close(pb)
  print("EM algorithm finished.")

  return(list(
    "mu_r" = mu_r, "Sigma_r" = Sigma_r, "beta_y" = beta_y, "sigma_y_sq" = sigma_y_sq, "beta_tte" = beta_tte, "sigma_tte_sq" = sigma_tte_sq, "pi_c" = pi_c,
    "mu_r_cure" = mu_r_cure, "Sigma_r_cure" = Sigma_r_cure, "beta_cure" = beta_cure_aug, "sigma_y_cure_sq" = sigma_y_cure_sq,
    "E_Delta_cen" = E_Delta_cen, "E_t_cen" = E_t_cen, "E_t_sq_cen" = E_t_sq_cen,
    "E_r_obs" = E_r_obs, "E_r_rT_obs" = E_r_rT_obs, "E_r_cen" = E_r_cen, "E_r_rT_cen" = E_r_rT_cen, "E_r_rT_mu_obs" = E_r_rT_mu_obs, "E_r_rT_mu_cen" = E_r_rT_mu_cen,
    "Estep_2_obs" = Estep_2_obs, "Estep_2_cen" = Estep_2_cen,
    "E_g0_t_cen" = E_g0_t_cen, "E_g1_t_cen" = E_g1_t_cen, "E_g2_t_cen" = E_g2_t_cen
  ))
}