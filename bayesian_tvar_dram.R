#tar_prior: list with distr and input as elements
bayesian_tvar <- function(data, P, tar_variable, max_d, tar_prior, tar_scale, tar_transform = NULL, 
                          reps = 10000, burn = reps / 2, stability = T,  max_attempts = 1e4, deltas = NULL, 
                          lambda = 0.1, tau = 10 * lambda, eps = 1e-4,  sum_of_coeff = T, 
                          train_sample = NULL, quiet = F, forecast_horizon = 24, irf_settings = NULL, 
                          dram_settings = list(start = 100, adapt = 50, adapt_scale = 2.4^2, 
                                               drscale = 2, delta = 1e-5)) {
  #ar_coeff = F,  #taken out from function call
  
  #tar scale is for the random walk MH; while the variance of the prior is in tar_prior
  #tar scale and tar prior provide variances -> sqrt to sd is done !!!!!!!!
  out_call <- match.call()
  if (reps <= burn)
    stop("Argument 'reps' has to be bigger than argument 'burn'.")
  if (is.character(tar_variable)) {
      tar_variable_character <- tar_variable
      tar_variable <- which(colnames(data) == tar_variable) 
      if (length(tar_variable) == 0)
        stop("Variable `", tar_variable_character, "' not found.")
  }
  N <- ncol(data)
  size <- N * P #for Bc_template and df for iW, also for computing n_crit and B_sample
  zero_vec <- rep(0, N) #for evaluating likelihoods with residuals -> expected value is zero vector
  #====================================
  # stdev <- F #for dummy prior
  # ar_p <- F #for dummy prior
  # ar_coeff <- T
  #===================================
  tar_prior$distr <- as.character(tar_prior$distr)
  T_complete <- nrow(data)
  d <- 1:max_d #all possible d values; current d will be saved in d_sample
  #Transform
  tardata <- as.matrix(data[, tar_variable, drop = F]) #later to be transformed maybe
  if (is.null(tar_transform))
    tar_transform <- list(fun = "ident", inp = NA)
  tardata <- tar_fn_cpp(tardata, tar_transform) #transform
  T_tardata <- nrow(tardata)
  tar_scale_adpat <- tar_scale / dram_settings$drscale #tar scale is still the VARIANCE! and so is tar_scale_adpat; done differently in DRAM: chol(), than div by 2(3)
  tar_scale <- sqrt(tar_scale) #for rnorm adjusted to sd
  tar_scale_adpat <- sqrt(tar_scale_adpat) #for rnorm adjusted to sd
  dtar <- match.fun(paste0("d", tar_prior$distr)) #density of tar variable
  if (tar_prior$distr == "norm") {
    tar_inp1 <- mean(tardata, na.rm = T) #tardata is not lagged yet!
    tar_inp2 <- sqrt(as.integer(tar_prior$inp)) #variance -> sd
  } 
  if (tar_prior$distr == "unif") { #STILL TO TEST
    if (is.null(tar_prior$inp)) {
      tar_inp1 <- min(tardata)
      tar_inp2 <- max(tardata)
    } else {
      tar_inp1 <- tar_prior$inp[1]
      tar_inp2 <- tar_prior$inp[2]
    }
  }
  
  #lag tar independently because of possible tranformation
  ZZ <- embed(tardata, max_d + 1) #creating threshold values
  ZZ <- ZZ[complete.cases(ZZ), -1, drop = F] #getting rid of current value
  data_embed <- embed(as.matrix(data), P + 1)
  diff_rows <- nrow(data_embed) - nrow(ZZ)
  #adjusting number of rows for ZZ or data_embed
  if (diff_rows < 0) 
    ZZ <- ZZ[-c(1:abs(diff_rows)), ]
  if (diff_rows > 0) 
    data_embed <- data_embed[-c(1:diff_rows), ]
  #delete train samples: STILL TO TEST; BUT NOT DONE IN ANY PAPERS -> COMMENTED OUT
  #------------------------------------
  # if (train_sample != T_complete) {
  #   ZZ <- ZZ[-c(1:train_sample), ] 
  #   data_embed <- data_embed[-c(1:train_sample), ] 
  # }
  #------------------------------------
  T_embed <- nrow(data_embed)
  YY <- data_embed[, 1:N]
  XX <- cbind(data_embed[, -c(1:N)], 1)
  n_crit <- size + 1 #min number of observations in each regime
  if (!stability)
    max_attempts <- 0 #no eigenvalues will be computed in B_sample, i.e. no stability test
  #dummy Prior
  if (is.null(train_sample))
    train_sample <- T_complete
  # dummy_p <- create_dummy_prior(data[1:train_sample, ], P, lambda, tau, eps, ar_coeff, ar_p,
  #                               stdev, sum_of_coeff)
  #DELTAS ARE ALL EQUAL TO ZERO BY DEFAULT!!!! 
  if (is.null(deltas))
    deltas <- rep(0, N)
  #NEXT ONE IS STANDARD
  dummy_p <- create_dummy_prior_man(data[1:train_sample, ], P = P, deltas = deltas, lambda = lambda,
                                     eps = eps, sum_of_coeff = sum_of_coeff, tau = tau)
  # mus <- colMeans(YY)
  # dummy_p <- create_dummy_prior_man2(data[1:train_sample, ], P = P, deltas = deltas, lambda = lambda, 
  #                                    eps = eps, sum_of_coeff = sum_of_coeff, tau = tau, mus = mus)
  YD <- dummy_p$Yd
  XD <- dummy_p$Xd
  
  #starting values
  d_sample <- sample(d, 1)
  curr_ZZ <- ZZ[, d_sample]
  #starting tar value
  curr_ZZ_sorted <- sort(curr_ZZ)
  curr_ZZ_sorted <- curr_ZZ_sorted[-c(1:(n_crit - 1), (T_embed - n_crit + 1):T_embed)] #guarantees number of obs sufficient
  tar_value <- sample(curr_ZZ_sorted, 1)
  reg1 <- curr_ZZ <= tar_value
  reg2 <- !reg1
  Y1 <- rbind(YY[reg1, ], YD); X1 <- rbind(XX[reg1, ], XD)
  Y2 <- rbind(YY[reg2, ], YD); X2 <- rbind(XX[reg2, ], XD)
  fit_reg1 <- .lm.fit(X1, Y1)
  fit_reg2 <- .lm.fit(X2, Y2)
  B1_sample <- c(fit_reg1$coefficients) #just in case we do not get stable results for B1, B2 in first Gibbs iteration
  B1_sample_mat <-  t(fit_reg1$coefficients)
  B2_sample <- c(fit_reg2$coefficients)
  B2_sample_mat <- t(fit_reg2$coefficients)
  Sigma1_lm <- crossprod(fit_reg1$residuals) #/ nrow(Y1)
  Sigma2_lm <- crossprod(fit_reg2$residuals) #  / nrow(Y2)
  repeat { #does that work?
    # Sigma1_sample <- riwish_cpp(nrow(Y1) + 2 - (size + 1), Sigma1_lm)  #random starting value; nrow(Y1) is after dummies appended -> correct
    # Sigma2_sample <- riwish_cpp(nrow(Y2) + 2 - (size + 1), Sigma2_lm)  
    Sigma1_sample <- riwish_cpp(nrow(Y1), Sigma1_lm)  #random starting value; nrow(Y1) is after dummies appended -> correct
    Sigma2_sample <- riwish_cpp(nrow(Y2), Sigma2_lm)  
    # Sigma1_sample <- diag(N); Sigma2_sample <- diag(N)
    if (min(eigen(Sigma1_sample, only.values = TRUE)$values, 
            eigen(Sigma2_sample, only.values = TRUE)$values) > 0)
      break
  }
  out_start <- list(Sigma1_start = Sigma1_sample, Sigma2_start = Sigma2_sample,  #starting values for return
                    start_tar = tar_value, start_d = d_sample)
  
  if (!quiet) 
    print(paste0("Starting TVAR estimation with ", reps, " replications."), quote = F)
  out_beta1 <- matrix(NA_real_, reps - burn, length(fit_reg1$coefficients))
  out_beta2 <- matrix(NA_real_, reps - burn, length(fit_reg2$coefficients))
  out_sigma1 <- array(NA_real_, dim = c(reps - burn, dim(Sigma1_sample)))
  out_sigma2 <- array(NA_real_, dim = c(reps - burn, dim(Sigma2_sample)))
  out_resid <- array(NA_real_, dim = c(reps - burn, dim(YY)))
  if (!is.null(forecast_horizon)) {
    forecast_start_period <- T_embed - P + 1  #first value to be used in forecast (as lag P)
    out_yhat <- array(NA_real_, dim = c(reps - burn,  P + forecast_horizon, N))
    start_forecast <- YY[forecast_start_period:T_embed, ]
  }
  out_tar <- rep(NA_real_, reps) #has to be reps, not reps - burn for DRAM procedure; later size is reduced before returned
  out_delay <- rep(NA_real_, reps - burn)
  out_post <- rep(NA_real_, reps - burn)
  if (!is.null(irf_settings)) {
    shocked_variable <- irf_settings$shocked_variable
    shock_size <- irf_settings$shock_size
    irf_horizon <- irf_settings$horizon
    restrict <- irf_settings$restrict
    type <- irf_settings$type
    out_ir1 <- array(NA_real_, dim = c(reps - burn, irf_horizon, N))
    out_ir2 <- array(NA_real_, dim = c(reps - burn, irf_horizon, N))
  }
  n_accept <- 0
  adjust_runs <- min(reps * .75, burn * .95) #when to stop the AM adjustment
  Bc_template <- rbind(matrix(0, N, size), cbind(diag(N * (P - 1)), matrix(0, N * (P - 1), N)))
  for (iter in 1:reps) {
    #sample B1 and Sigma1
    reg1 <- curr_ZZ <= tar_value
    reg2 <- !reg1
    T_reg1 <- sum(reg1) #for deleting the residuals from prior when finding new thresh value (see eval_post)
    T_reg2 <- sum(reg2)
    Y1 <- rbind(YY[reg1, ], YD)
    X1 <- rbind(XX[reg1, ], XD)
    B1_star <- c(.lm.fit(X1, Y1)$coefficients)
    xpx1_inv <- xpx_inv_cpp(X1)
    B_return <- sample_B_fast(B1_star, Sigma1_sample, xpx1_inv, Bc_template, size, N, P, max_attempts)
    chk1 <- B_return$chk #stable
    if (chk1 || !stability) {
      B1_sample <- B_return$B_sample
      B1_sample_mat <- B_return$B_mat
    }
    resid1 <- resids_cpp(Y1, X1, B1_sample_mat) #all resids, including from dummies for next step
    # Sigma1_sample <- riwish_cpp(nrow(Y1) + 2 - (size + 1), crossprod(resid1))
    Sigma1_sample <- riwish_cpp(nrow(Y1), crossprod(resid1))
    
    #sample B2 and Sigma2
    Y2 <- rbind(YY[reg2, ], YD)
    X2 <- rbind(XX[reg2, ], XD)
    B2_star <- c(.lm.fit(X2, Y2)$coefficients)
    xpx2_inv <- xpx_inv_cpp(X2)
    B_return <- sample_B_fast(B2_star, Sigma2_sample, xpx2_inv, Bc_template, size, N, P, max_attempts)
    chk2 <- B_return$chk #stable
    if (chk2 || !stability) {
      B2_sample <- B_return$B_sample
      B2_sample_mat <- B_return$B_mat
    }
    resid2 <- resids_cpp(Y2, X2, B2_sample_mat)
    # Sigma2_sample <- riwish_cpp(nrow(Y2) + 2 - (size + 1), crossprod(resid2))
    Sigma2_sample <- riwish_cpp(nrow(Y2), crossprod(resid2))
 
    #sample tar: MH step
    tar_value_star <- rnorm(1, tar_value, tar_scale) #sample new tar value )
    post_old <- dtar(tar_value, tar_inp1, tar_inp2, log = T) + 
      lik_cpp(resid1[1:T_reg1, ], zero_vec, Sigma1_sample, loglik = T) + 
      lik_cpp(resid2[1:T_reg2, ], zero_vec, Sigma2_sample, loglik = T)
    reg1 <- curr_ZZ <= tar_value_star
    reg2 <- !reg1
    Y1_new <- YY[reg1,, drop = F]; X1_new <- XX[reg1,, drop = F]
    Y2_new <- YY[reg2,, drop = F]; X2_new <- XX[reg2,, drop = F]
    if (min(nrow(Y1_new), nrow(Y2_new)) >= n_crit) { #also done in matlab by setting post to -Inf if nobs < n_crit, so never accepted
      resid1 <- resids_cpp(Y1_new, X1_new, B1_sample_mat)
      resid2 <- resids_cpp(Y2_new, X2_new, B2_sample_mat)
      post_new <- dtar(tar_value_star, tar_inp1, tar_inp2, log = T) + #don't have to delete "dummy resids" because not included here
        lik_cpp(resid1, zero_vec, Sigma1_sample, loglik = T) + 
        lik_cpp(resid2, zero_vec, Sigma2_sample, loglik = T)
      alpha_12 <- min(1, exp(post_new - post_old)) #12...moving from 1 to 2
      if (alpha_12 > runif(1)) {
        tar_value <- tar_value_star
        n_accept <- n_accept + 1
      } else {#if (iter < adjust_runs) {#we didn't accept the first attempt
        tar_value_star_2 <- rnorm(1, tar_value, tar_scale_adpat) #sample new tar value with smaller scale
        reg1 <- curr_ZZ <= tar_value_star_2
        reg2 <- !reg1
        Y1_new <- YY[reg1,, drop = F]; X1_new <- XX[reg1,, drop = F]
        Y2_new <- YY[reg2,, drop = F]; X2_new <- XX[reg2,, drop = F]
        if (min(nrow(Y1_new), nrow(Y2_new)) >= n_crit) {
          resid1 <- resids_cpp(Y1_new, X1_new, B1_sample_mat)
          resid2 <- resids_cpp(Y2_new, X2_new, B2_sample_mat)
          post_new_2 <- dtar(tar_value_star_2, tar_inp1, tar_inp2, log = T) +
            lik_cpp(resid1, zero_vec, Sigma1_sample, loglik = T) +
            lik_cpp(resid2, zero_vec, Sigma2_sample, loglik = T)
          alpha_32 <- min(1, exp(post_new - post_new_2)) #moving from 3 to 2
          q_ratio <- dnorm(tar_value_star, tar_value_star_2, tar_scale, log = T) -
            dnorm(tar_value_star, tar_value, tar_scale, log = T) #moving from star_2 to star compared to moving from current to star
          ll_ratio <- post_new_2 - post_old
          alpha_13 <- exp(ll_ratio + q_ratio) * (1 - alpha_32) / (1 - alpha_12)
          if (alpha_13 > runif(1)) {
            tar_value <- tar_value_star_2
            n_accept <- n_accept + 1
          }
        }
      }
    }
    a_rate <- n_accept / iter
    out_tar[iter] <- tar_value
    
    #cat(tar_scale, " ", a_rate, "\n")
    #adjust tar scale
    if (iter < adjust_runs && iter >= dram_settings$start && iter %% dram_settings$adapt == 0) {
      tar_scale_new <- var(out_tar[1:iter]) + dram_settings$delta #find variance, but we need sd!
      #tar_scale_new <- var(out_tar[1:iter]) * dram_settings$adapt_scale
      if (tar_scale_new != 0) {
        tar_scale <- tar_scale_new * dram_settings$adapt_scale #variance
        tar_scale_adpat <- tar_scale / dram_settings$drscale #variance
        tar_scale <- sqrt(tar_scale) #make sd
        tar_scale_adpat <- sqrt(tar_scale_adpat) #make sd
      }
    }
    
    #sample delay
    probs <- eval_delay_thresh_cpp(ZZ, tar_value, YY, XX, B1_sample_mat, B2_sample_mat, 
                                   Sigma1_sample, Sigma2_sample)
    d_sample <- sample(d, 1, prob = probs)
    curr_ZZ <- ZZ[, d_sample]
    
    # if (!quiet && iter == (burn + 1))
    #   print("Burn in phase done.", quote = F)
    if (iter > burn) {
      out_beta1[iter - burn, ] <- B1_sample
      out_beta2[iter - burn, ] <- B2_sample
      out_sigma1[iter - burn,, ] <- Sigma1_sample
      out_sigma2[iter - burn,, ] <- Sigma2_sample
      out_delay[iter - burn] <- d_sample
      #saving post value
      reg1 <- curr_ZZ <= tar_value
      reg2 <- !reg1
      Y1_new <- YY[reg1, ]; X1_new <- XX[reg1, ]
      Y2_new <- YY[reg2, ]; X2_new <- XX[reg2, ]
      resid1 <- resids_cpp(Y1_new, X1_new, B1_sample_mat)
      resid2 <- resids_cpp(Y2_new, X2_new, B2_sample_mat)
      out_resid[iter - burn, which(reg1), ] <- resid1 #correct order now
      out_resid[iter - burn, which(reg2), ] <- resid2
      # resid1 <- Y1_new - X1_new %*% t(B1_sample_mat)
      # resid2 <- Y2_new - X2_new %*% t(B2_sample_mat)
      out_post[iter - burn] <-  dtar(tar_value, tar_inp1, tar_inp2, log = T) + 
        lik_cpp(resid1, zero_vec, Sigma1_sample, loglik = T) + 
        lik_cpp(resid2, zero_vec, Sigma2_sample, loglik = T)
      if (((chk1 && chk2) || !stability) && !is.null(forecast_horizon)) {
        out_yhat[iter - burn,, ] <- tvar_forecast_cpp(start_forecast, tar_variable, P, N, 
                                                      forecast_horizon, d_sample,T_tardata, 
                                                      tar_value, tardata,B1_sample_mat,
                                                      B2_sample_mat, Sigma1_sample,
                                                      Sigma2_sample, tar_transform, data)
      }
      if (((chk1 && chk2) || !stability) && !is.null(irf_settings)) {
        out_ir1[iter - burn,, ] <- irf(type, B1_sample_mat, Sigma1_sample, shocked_variable,
                                       shock_size, irf_horizon, restrict) 
        out_ir2[iter - burn,, ] <- irf(type, B2_sample_mat, Sigma2_sample, shocked_variable,
                                       shock_size, irf_horizon, restrict) 
      }
    }
    if (!quiet && iter %% 1000 == 0)
      print(paste0("Replication ", iter, " of ", reps, ". Acceptance Ratio = ", 
                   round(a_rate, 5), ". Sqrt tarscale = ", round(tar_scale, 5), "."), quote = F)
  }
  #CHANGED out_yhat!!!!!!
  ret_list <- list(out_beta1 = out_beta1, out_beta2 = out_beta2, out_sigma1 = out_sigma1,
                   out_sigma2 = out_sigma2, out_yhat = NULL, out_tar = out_tar[-c(1:burn)], 
                   out_delay = out_delay, out_post = out_post, out_resid = out_resid)
  if (!is.null(irf_settings)) {
    ret_list$out_ir1 <- out_ir1
    ret_list$out_ir2 <- out_ir2
  }
  if (!is.null(forecast_horizon))
    ret_list$out_yhat <- out_yhat
  ret_list$acceptance_rate <- a_rate
  ret_list$tar_scale <- tar_scale^2  #CHANGED TO ^2 SO WE RETURN VARIANCE AGAIN (input is variance and then transformed to sd)
  ret_list$starting_values <- out_start
  ret_list$model_specific <- list(ZZ = ZZ, data_embed = data_embed, dummy_p = dummy_p)
  ret_list$out_call <- out_call
  ret_list$dataset <- data
  #ret_list$all_probs <- all_probs
  return(ret_list)
}



