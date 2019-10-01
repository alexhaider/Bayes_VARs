bayesian_var <- function(data, P, reps = 10000, burn = 5000, stability = T, max_attempts = 1e4,
                         deltas = NULL, lambda = 0.1, tau = 10 * lambda, eps = 1e-4, 
                         sum_of_coeff = T, train_sample = NULL, quiet = F, forecast_horizon = 24, 
                         irf_settings = NULL) {
  #ar_coeff = F,  #taken out from function calls
  out_call <- match.call()
  if (reps <= burn) 
    stop("Argument 'reps' has to be bigger than argument 'burn'.")
  N <- ncol(data)
  size <- N * P #for Bc_template, also used for df for iW distribution: df = T_star + 2 - (size + 1)
  T_complete <- nrow(data) #needed if no value for training sample provided (standard approach)
  ############
  # stdev <- F #for dummy prior
  # ar_p <- F
  # ar_coeff <- T
  ############
  if (!stability)
    max_attempts <- 0 #no eigenvalues computed in sample_B -> no stability test
  if (is.null(train_sample))
    train_sample <- T_complete
  # dummy_p <- create_dummy_prior(data[14:train_sample, ], P, lambda, tau, eps, ar_coeff, F, sum_of_coeff)
  # dummy_p <- create_dummy_prior(data[1:train_sample, ], P, lambda, tau, eps, ar_coeff, ar_p,
  #                              stdev, sum_of_coeff)
  if (is.null(deltas))
    deltas <- rep(0, N) #assumption of transfomred/stationary data if no deltas provided
  #NEXT ONE IS STANDARD
  dummy_p <- create_dummy_prior_man(data[1:train_sample, ], P = P, deltas = deltas, lambda = lambda,
                                    eps = eps, sum_of_coeff = sum_of_coeff, tau = tau)
  # mus <- colMeans(data[-c(1:P), ])
  # dummy_p <- create_dummy_prior_man(data[1:train_sample, ], P = P, deltas = deltas, lambda = lambda, 
  #                                   eps = eps, sum_of_coeff = sum_of_coeff, tau = tau, mus = mus)
  YD <- dummy_p$Yd
  XD <- dummy_p$Xd
  data_embed <- embed(as.matrix(data), P + 1)
  T_embed <- nrow(data_embed) #effective number of observation afer lagging
  #REMOVE TRAININGS SAMPLE: STILL TO TEST; BUT NOT DONE IN ANY PAPERS -> COMMENTED OUT
  # if (train_sample != T_complete) {
  #   YY <- data_embed[-c(1:train_sample), 1:N] #last lag in X is first obs after training sample
  #   XX <- cbind(data_embed[-c(1:train_sample), -c(1:N)], 1)
  # } else {
  #   YY <- data_embed[, 1:N]
  #   XX <- cbind(data_embed[, -c(1:N)], 1)
  # }
  YY <- data_embed[, 1:N] #new because remove training sample has been commented out
  XX <- cbind(data_embed[, -c(1:N)], 1) #new because remove training sample has been commented out
  Y_star <- rbind(YY, YD)
  X_star <- rbind(XX, XD)
  #OLS estimate and other prep
  T_star <- nrow(Y_star) #nobs with dummy obs: Tstar = T_embed + T_d -> T_star for df of iW distribution
  xpx_inv <- xpx_inv_cpp(X_star)
  B_star_mat <- .lm.fit(X_star, Y_star)$coefficients
  B_star <- c(B_star_mat)
  #starting value
  S_star <- resids_cpp(Y_star, X_star, t(B_star_mat), crossprod = TRUE) #works as expected!
  # Sigma_sample <- riwish_cpp(T_star + 2 - (size + 1), S_star) #T_star = T_d + T
  Sigma_sample <- riwish_cpp(T_star, S_star)
  # if (!quiet) 
  #   pb <- txtProgressBar(min = 0, max = reps, style = 3)
  if (!quiet) 
    print(paste0("Starting VAR estimation with ", reps, " replications."), quote = F)
  out_beta <- matrix(NA_real_, reps - burn, length(B_star))
  out_sigma <- array(NA_real_, dim = c(reps - burn, dim(S_star)))
  out_resid <- array(NA_real_, dim = c(reps - burn, dim(YY)))
  if (!is.null(forecast_horizon)) {
    forecast_start_period <- T_embed - P + 1 
    forecast_start <- YY[forecast_start_period:T_embed, ]
    out_yhat <- array(NA_real_, dim = c(reps - burn,  P + forecast_horizon, N))
  }
  if (!is.null(irf_settings)) {
    shocked_variable <- irf_settings$shocked_variable
    shock_size <- irf_settings$shock_size
    irf_horizon <- irf_settings$horizon
    restrict <- irf_settings$restrict
    type <- irf_settings$type
    out_ir <- array(NA_real_, dim = c(reps - burn, irf_horizon, N))
  }
  #Companion Matrix template
  Bc_template <- rbind(matrix(0, N, size), cbind(diag(N * (P - 1)), matrix(0, N * (P - 1), N)))
  for (iter in 1:reps) {
    #cat(iter, "\n")
    B_return <- sample_B_fast(B_star, Sigma_sample, xpx_inv, Bc_template, size, N, P, max_attempts)
    chk <- B_return$chk #stable
    if (chk || !stability) {
      B_sample <- B_return$B_sample
      B_sample_mat <- B_return$B_mat
    }
    S_star <- resids_cpp(Y_star, X_star, B_sample_mat, crossprod = TRUE)
    # Sigma_sample <- riwish_cpp(T_star + 2 - (size + 1), S_star)
    Sigma_sample <- riwish_cpp(T_star, S_star)
    if (iter > burn) {
      out_beta[iter - burn, ] <- B_sample
      out_sigma[iter - burn,, ] <- Sigma_sample
      out_resid[iter - burn,, ] <- resids_cpp(Y_star, X_star, B_sample_mat, 
                                              crossprod = FALSE)[1:T_embed, ] #1:T_embed to get rid of dummy obs resids
      if ((chk || !stability) && !is.null(forecast_horizon)) {
        out_yhat[iter - burn,, ] <- lin_forecast_cpp(forecast_start, P, N, forecast_horizon, 
                                                     B_sample_mat, Sigma_sample)
        if (!is.null(irf_settings)) 
          out_ir[iter - burn,, ] <- irf(type, B_sample_mat, Sigma_sample, shocked_variable,
                                        shock_size, irf_horizon, restrict) 
      }
    }
    if (!quiet && iter %% 1000 == 0)
      print(paste0("Replication ", iter, " of ", reps, "."), quote = F)
    # if (!quiet)
    #   setTxtProgressBar(pb, iter)
  }
  est_res <- list(out_beta = out_beta, out_sigma = out_sigma, out_resid = out_resid)
  if (!is.null(forecast_horizon))
    est_res$out_yhat <- out_yhat
  if (!is.null(irf_settings)) 
    est_res$out_ir <- out_ir
  est_res$model_specific <- list(data_embed = data_embed, P = P)
  est_res$out_call <- out_call
  est_res$dataset <- data
  return(est_res)
}

