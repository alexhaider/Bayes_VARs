#go and Tg0 for stoch vol error term, V0 = variance for A matrix
bayesian_var_svol <- function(data, P, sigbar = 10, g0 = 1e-4, Tg0 = 1, V0 = 1000, 
                              reps = 10000, burn = reps / 2, stability = T, max_attempts = 1e3, 
                              deltas = NULL, lambda = 0.1, tau = 10 * lambda, eps = 1e-4, 
                              sum_of_coeff = T, train_sample = 40, quiet = F, 
                              forecast_horizon = 24, irf_settings = NULL) {
  #ar_coeff = F  #taken out from function call
  #deterministic starting values for hlast and res
  out_call <- match.call()
  if (reps <= burn)
    stop("Argument 'reps' has to be bigger than argument 'burn'.")
  #Stability Testing?
  if (!stability)
    max_attempts <- 0 #should work given construction of CK algo
  N <- ncol(data)
  #size <- N * P
  #diag_mat <- diag(N)
  
  #stdev <- TRUE #for dummy prior
  #ar_p <- TRUE
  #Prior for Beta as starting value for CK ==
  # dummy_p <- create_dummy_prior(data[-c(1:P), ], P, lambda, tau, eps, ar_coeff, ar_p,
  #                               stdev, sum_of_coeff)
  
  #new Dummy Prior: second version again
  # mus <- colMeans(data[-c(1:P), ])
  # dummy_p <- create_dummy_prior_man2(data[1:train_sample, ], P = P, deltas = deltas, lambda = lambda, 
  #                                    eps = eps, sum_of_coeff = sum_of_coeff, tau = tau, mus = mus)
  if (is.null(deltas)) 
    deltas <- rep(0, N)
  #data instead of data[1:train_sample, ] is okay in prior because we set train_sample = nobs anyway...
  dummy_p <- create_dummy_prior_man(data, P = P, deltas = deltas, lambda = lambda, eps = eps, 
                                    sum_of_coeff = sum_of_coeff, tau = tau) #used for starting vals for Kalman filter
 
  #YD <- dummy_p$Yd[-c((P * N + 1):(P * N + 4)), ] #remains in there?
  #XD <- dummy_p$Xd[-c((P * N + 1):(P * N + 4)), ]
  lm_fit <- .lm.fit(dummy_p$Xd, dummy_p$Yd)
  b0_ck <- matrix(lm_fit$coefficients, nrow = 1) #start value for Carter Kohn
  #size_b <- length(b0_ck) #NEW
  Sig <- diag(dummy_p$sig^2)
  s0 <- kronecker(Sig, xpx_inv_cpp(dummy_p$Xd)) #start value for Carter Kohn #different now
  
  #Prepare data ==
  data_embed <- embed(data, P + 1)
  #nobs <- nrow(data_embed) #eff_nobs below #nobs has to be evaluated due to error in Matlab Code when forming prior for hlast 
  YY <- data_embed[, 1:N]
  XX <- cbind(data_embed[, -c(1:N)], 1)
  YY <- YY[-c(1:(train_sample - P)), ] #GUARANTEES NO OVERLAP WITH Y0 BELOW; SHOULDN'T IT BE TRAIN_SAMPLE - P
  XX <- XX[-c(1:(train_sample - P)), ] #GUARANTEES NO OVERLAP WITH Y0 BELOW
  T_embed <- nrow(YY)
  lm_xy <- .lm.fit(XX, YY) #starting values for residuals and beta (although beta not really used)
  beta <- matrix(c(lm_xy$coefficients), nrow = 1) #initialization of beta; not needed in algorithm but still...
  res <- lm_xy$residuals #residual starting values
  
  #Prior, starting values for hlast and g ==
  presample <- embed(data, 2)
  Y0 <- presample[1:(train_sample - 1), 1:N] #in Matalb: train_sample - 1
  X0 <- cbind(presample[1:(train_sample - 1), -c(1:N)], 1) #in Matalb: train_sample - 1
  lm_fit0 <- .lm.fit(X0, Y0)
  resid0 <- lm_fit0$residuals
  Sigma0 <- crossprod(resid0) / (train_sample - 1)  #should be T0!? used for mubar prior only #was nobs
  mubar <- log(diag(Sigma0)) #prior for h_0 (mean), variance is set provided by sigbar (function call)
  hlast <- diff(YY)^2
  #!!!! sd(y) / 1000? instead like Kastner??
  hlast <- rbind(hlast[1, ], hlast, hlast[nrow(hlast), ]) + 0.0001 
  #all.equal(data[-1, ], rbind(Y0, YY)) #first obs has to be deleted because of lagging, no overlap, no info lost
  
  #out_start <- list(Sigma1_start = Sigma1_sample, Sigma2_start = Sigma2_sample,  #starting values for return
  #                  start_tar = tar_value, start_d = d_sample)
  
  if (!quiet) 
    print(paste0("Starting BVAR_SVOL estimation with ", reps, " replications."), quote = F)
  #Q <- matrix(0, (N * (N * P + 1)), (N * (N * P + 1)))
  #mu_beta <- rep(0, size_b) 
  #Fmat_beta <- diag(size_b)
  write_iter <- 1
  out_A <- matrix(NA_real_, reps - burn, N * (N - 1) / 2)
  out_g <- matrix(NA_real_, reps - burn, N)
  out_B <- matrix(NA_real_, reps - burn, N^2 * P + N)
  out_h <- array(NA_real_, dim = c(reps - burn, T_embed, N))
  out_resid <- array(NA_real_, dim = c(reps - burn, dim(YY)))
  if (!is.null(forecast_horizon)) {
    forecast_start_period <- T_embed - P + 1  #first value to be used in forecast (as lag P)
    out_yhat <- array(NA_real_, dim = c(reps - burn,  P + forecast_horizon, N))
    start_forecast <- YY[forecast_start_period:T_embed, ]
  }
  # if (!is.null(irf_settings)) {
  #   shocked_variable <- irf_settings$shocked_variable
  #   shock_size <- irf_settings$shock_size
  #   irf_horizon <- irf_settings$horizon
  #   restrict <- irf_settings$restrict
  #   type <- irf_settings$type
  #   out_ir1 <- array(NA_real_, dim = c(reps - burn, irf_horizon, N))
  #   out_ir2 <- array(NA_real_, dim = c(reps - burn, irf_horizon, N))
  # }
  Bc_template <- rbind(matrix(0, N, N * P), cbind(diag(N * (P - 1)), matrix(0, N * (P - 1), N)))
  H_cube <- comp_Hmat(XX, N) # for Kalman filter: eye %*% X 
  for (iter in 1:reps) {
    #sample g
    g <- sample_g(hlast, Tg0, g0)
    if (any(is.na(g)))
      stop("Variance 'g' is negative for at least one element: ", g)
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
    hlast <- svol_mh(hlast, ortho_res, g, mubar, sigbar)
    
    if (iter > burn) {
      out_A[write_iter, ] <- a_complete$a_draw
      out_g[write_iter, ] <- g
      out_B[write_iter, ] <- beta
      out_h[write_iter,, ] <- hlast[-1, ] #delete starting value
      if (!is.null(forecast_horizon)) 
        out_yhat[write_iter,, ] <- svol_ca_forecast_cpp(start_forecast, P, N, forecast_horizon,
                                                        ck_draw$B_mat, hlast, g, A)
      out_resid[write_iter,, ] <- res
      write_iter <- write_iter + 1
    }
    
    if (!quiet && iter %% 1000 == 0)
      print(paste0("Replication ", iter, " of ", reps, "."), quote = F)
  }
  ret_list <- list(out_B = out_B, out_A = out_A, out_h = out_h, out_g = out_g, out_resid = out_resid)
  # if (!is.null(irf_settings)) {
  #   ret_list$out_ir1 <- out_ir1
  #   ret_list$out_ir2 <- out_ir2
  # }
  if (!is.null(forecast_horizon))
    ret_list$out_yhat <- out_yhat
  #ret_list$starting_values <- out_start
  ret_list$model_specific <- list(data_embed = data_embed, dummy_p = dummy_p)
  ret_list$out_call <- out_call
  ret_list$dataset <- data
  return(ret_list)
}