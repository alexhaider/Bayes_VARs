#go and Tg0 for stoch vol error term, V0 = variance for A matrix
bayesian_var_svol_ar <- function(data, P, sv_priors, h_start = NULL, V0 = 1000, reps = 10000, 
                                 burn = reps / 2, stability = T, max_attempts = 1e3, deltas = NULL, 
                                 lambda = 0.1, tau = 10 * lambda, eps = 1e-4, sum_of_coeff = T, 
                                 train_sample = 40, quiet = F, forecast_horizon = 24, 
                                 irf_settings = NULL) {
  #deterministic starting values for hlast and res
  out_call <- match.call()
  if (reps <= burn)
    stop("Argument 'reps' has to be bigger than argument 'burn'.")
  #Stability Testing?
  if (!stability)
    max_attempts <- 0 #should work given construction of CK algo
  N <- ncol(data)
  
  #Starting values for CK ==
  if (is.null(deltas)) 
    deltas <- rep(0, N)
  #data instead of data[1:train_sample, ] is okay in prior because we set train_sample = nobs anyway...
  dummy_p <- create_dummy_prior_man(data, P = P, deltas = deltas, lambda = lambda, eps = eps, 
                                    sum_of_coeff = sum_of_coeff, tau = tau) #used for starting vals for Kalman filter
  lm_fit <- .lm.fit(dummy_p$Xd, dummy_p$Yd)
  b0_ck <- matrix(lm_fit$coefficients, nrow = 1) #start value for Carter Kohn
  #size_b <- length(b0_ck) #NEW
  Sig <- diag(dummy_p$sig^2)
  s0 <- kronecker(Sig, xpx_inv_cpp(dummy_p$Xd)) #start value for Carter Kohn #different now
  
  #Prepare data and intialize residuals (and beta) ==
  data_embed <- embed(data, P + 1)
  #nobs <- nrow(data_embed) #eff_nobs below #nobs has to be evaluated due to error in Matlab Code when forming prior for hlast 
  YY <- data_embed[, 1:N]
  XX <- cbind(data_embed[, -c(1:N)], 1)
  T_embed <- nrow(YY)
  lm_xy <- .lm.fit(XX, YY) #starting values for residuals and beta (although beta not really used)
  beta <- matrix(c(lm_xy$coefficients), nrow = 1) #initialization of beta; not needed in algorithm but still...
  res <- lm_xy$residuals #residual starting values
  
  #Priors and Start for Stoch vol ==
  hlast <- diff(YY)^2
  hlast <- rbind(hlast[1, ], hlast, hlast[nrow(hlast), ]) + 0.0001
  if (is.null(h_start)) {
    h_start <- list()
    for (i in 1:N) 
      h_start[[i]] <- list(para = c(mu = 0, phi = .9, sigma = .1), latent = hlast[-1, i]) #starting
      # h_start[[i]] <- list(para = c(mu = 0, phi = .9, sigma = .1), latent = rep(0, T_embed)) #starting
  }
  hlast_all <- h_start #full info
  # hlast <- matrix(0, nrow = T_embed + 1, ncol = ncol(data)) #+1 to recycle code from svol with starting value for h0
  priormu <- sv_priors$priormu
  priorphi <- sv_priors$priorphi
  priorsigma <- sv_priors$priorsigma
  
  
  #Prior, starting values for hlast and g ==
  # presample <- embed(data, 2)
  # Y0 <- presample[1:(train_sample - 1), 1:N] #in Matalb: train_sample - 1
  # X0 <- cbind(presample[1:(train_sample - 1), -c(1:N)], 1) #in Matalb: train_sample - 1
  # lm_fit0 <- .lm.fit(X0, Y0)
  # resid0 <- lm_fit0$residuals
  # hlast <- diff(YY)^2
  # hlast <- rbind(hlast[1, ], hlast, hlast[nrow(hlast), ]) + 0.0001 
  
  if (!quiet) 
    print(paste0("Starting BVAR_SVOL_AR estimation with ", reps, " replications."), quote = F)
  write_iter <- 1
  out_A <- matrix(NA_real_, reps - burn, N * (N - 1) / 2)
  out_mu <- matrix(NA_real_, reps - burn, N)
  out_phi <- matrix(NA_real_, reps - burn, N)
  out_sig <- matrix(NA_real_, reps - burn, N)
  out_B <- matrix(NA_real_, reps - burn, N^2 * P + N)
  out_h <- array(NA_real_, dim = c(reps - burn, T_embed, N))
  out_resid <- array(NA_real_, dim = c(reps - burn, dim(YY)))
  if (!is.null(forecast_horizon)) {
    forecast_start_period <- T_embed - P + 1  #first value to be used in forecast (as lag P)
    out_yhat <- array(NA_real_, dim = c(reps - burn,  P + forecast_horizon, N))
    start_forecast <- YY[forecast_start_period:T_embed, ]
  }
  Bc_template <- rbind(matrix(0, N, N * P), cbind(diag(N * (P - 1)), matrix(0, N * (P - 1), N)))
  H_cube <- comp_Hmat(XX, N) # for Kalman filter: eye %*% X 
  for (iter in 1:reps) {
    #sample a
    a_complete <- sample_A(res, hlast, V0) #a_draw (vector) is used in out_A: easier saving vector
    A <- a_complete$Amat
    #Carter Kohn
    ck_draw <- carterkohn_ca_cb_fast2(YY, XX, A, hlast, b0_ck, s0, Bc_template, H_cube, max_attempts)
    # ck_draw <- carterkohn_ca_cb_fast(YY, XX, A, hlast, b0_ck, s0, Bc_template, max_attempts)
    #cat(ck_draw$attempts, " ")
    if (ck_draw$stability || stability == FALSE) 
      beta <- ck_draw$beta
    res <- ck_draw$errors
    ortho_res <- ck_draw$ortho_errors
    #sample h
    curr_mu <- curr_phi <- curr_sig <- rep(NA, N)
    for (svol_iter in 1:N) { #eps instead of res?
      h_curr <- hlast_all[[svol_iter]]
      curr_err <- ortho_res[, svol_iter]
      curr_err[which(curr_err == 0)] <- 1e-4 #no zeros for samplers
      hlast_all[[svol_iter]] <- svsample2(curr_err, startpara = para(h_curr), 
                                          startlatent = latent(h_curr), priormu = priormu, 
                                          priorphi = priorphi, priorsigma = priorsigma)
      hlast[-1, svol_iter] <- exp(latent(hlast_all[[svol_iter]])) #so that code for A,... works
      curr_para <- para(hlast_all[[svol_iter]]) 
      curr_mu[svol_iter] <- curr_para[1]
      curr_phi[svol_iter] <- curr_para[2]
      curr_sig[svol_iter] <- curr_para[3]
    }
    #saving the results
    if (iter > burn) {
      out_A[write_iter, ] <- a_complete$a_draw
      out_mu[write_iter, ] <- curr_mu
      out_phi[write_iter, ] <- curr_phi
      out_sig[write_iter, ] <- curr_sig
      out_B[write_iter, ] <- beta
      out_h[write_iter,, ] <- hlast[-1, ] #delete starting value
      if (!is.null(forecast_horizon))
        out_yhat[write_iter,, ] <- svol_ca_ar_forecast_cpp(start_forecast, P, N, forecast_horizon,
                                                           ck_draw$B_mat, hlast, curr_mu, curr_phi,
                                                           curr_sig, A)
      out_resid[write_iter,, ] <- res
      write_iter <- write_iter + 1
    }
    
    if (!quiet && iter %% 1000 == 0)
      print(paste0("Replication ", iter, " of ", reps, "."), quote = F)
  }
  ret_list <- list(out_B = out_B, out_A = out_A, out_h = out_h, out_mu = out_mu, out_phi = out_phi,
                   out_sig = out_sig, out_resid = out_resid)
  if (!is.null(forecast_horizon))
    ret_list$out_yhat <- out_yhat
  #ret_list$starting_values <- out_start
  ret_list$model_specific <- list(data_embed = data_embed, dummy_p = dummy_p)
  ret_list$out_call <- out_call
  ret_list$dataset <- data
  return(ret_list)
}