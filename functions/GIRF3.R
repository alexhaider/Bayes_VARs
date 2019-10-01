# est_res <- tvar_res1; shock_var <- 1; shock_size <- 1; type = "Triangular"
# horizon <- 20; n_histories <- 3; n_boot <- 4; quiet <- T
GIRF <- function(est_res, shock_var = 1, shock_size = 1, type = c("Chol", "Triangular", "Sign_restr"),
                 n_histories = 500, n_boot = 100, horizon = 20, quiet = FALSE) {
  #variables
  type <- match.arg(type)
  dataset <- eval(est_res$out_call$data) #complete dataset
  tar_variable <- eval(est_res$out_call$tar_variable)
  tar_transform <- eval(est_res$out_call$tar_transform)
  if (is.null(tar_transform) || tar_transform$fun == "ident")
    tar_transform <- list(fun = "ident", inp = 1)
  P <- eval(est_res$out_call$P)
  N <- ncol(dataset)
  YY <- est_res$model_specific$data_embed[, 1:N]
  XX <- cbind(est_res$model_specific$data_embed[, -c(1:N)], 1)
  ZZ <- est_res$model_specific$ZZ
  tardata <- dataset[, tar_variable]
  rows_to_del <- length(tardata) - (nrow(YY) + P + 1) #experimental: delete some rows because tardata is too long because of building MA,...
  if (rows_to_del > 0)
    tardata <- tardata[-c(1:rows_to_del)]
  if (rows_to_del < 0) #add some elements to make index work
    tardata <- c(rep(NA, abs(rows_to_del)), tardata)
  horizon <- horizon + 1
  n_regimes <- 2
  n_samples <- est_res$out_call$reps - est_res$out_call$burn
  if (is.character(shock_var)) {
    char_shock_var <- shock_var
    shock_var <- which(colnames(dataset) == shock_var)
    if (length(shock_var) == 0)
      stop("Error. Variable ", char_shock_var, " does not exist!")
  }
  
  #Loop over samples
  sample_iter <- 1
  final_girf_all_samples <- array(NA, dim = c(horizon, N, n_samples, n_regimes))
  for (sample_iter in 1:n_samples) {
    curr_B1 <- matrix(est_res$out_beta1[sample_iter, ], nrow = N, byrow = T)
    curr_B2 <- matrix(est_res$out_beta2[sample_iter, ], nrow = N, byrow = T)
    curr_B_array <- array(c(curr_B1, curr_B2), dim = c(nrow(curr_B1), ncol(curr_B1), n_regimes))
    curr_delay <- est_res$out_delay[sample_iter]
    max_delay_needed <- curr_delay + tar_transform$inp - 1 #transformation of TAR Variable -> more lags needed
    curr_tar <- est_res$out_tar[sample_iter]
    reg1 <- ZZ[, curr_delay] <= curr_tar; reg2 <- !reg1
    pos_in_ds <- list(reg1 = which(reg1), reg2 = which(reg2)) #saving postiiotns of obs in ds
    XX_1 <- XX[reg1,, drop = F]; XX_2 <- XX[reg2,, drop = F]
    resid_reg1 <- resids_cpp(YY[reg1, ], XX_1, curr_B1)
    resid_reg2 <- resids_cpp(YY[reg2, ], XX_2, curr_B2)
    
    #Prepare Sigma, histories and structural errors
    struc_errors <- histories <- list()
    struc_Sigma <- array(NA, dim = c(N, N, n_regimes))
    struc_Sigma[,, 1] <- struc_covmat_cpp(est_res$out_sigma1[sample_iter,, ], type) #lower triang
    struc_Sigma[,, 2] <- struc_covmat_cpp(est_res$out_sigma2[sample_iter,, ], type)
    struc_errors[[1]] <- struc_error_cpp(struc_Sigma[,, 1], resid_reg1) #inversion Sigma inside CPP function
    struc_errors[[2]] <- struc_error_cpp(struc_Sigma[,, 2], resid_reg2)
    histories[[1]] <- XX_1; histories[[2]] <- XX_2
    
    sim_res <- sim_sample_cpp(histories, struc_errors, pos_in_ds, horizon, N, P, max_delay_needed, 
                              curr_B_array, struc_Sigma, tardata, curr_tar, tar_variable, shock_var, 
                              shock_size, tar_transform, n_histories, n_boot)
    final_girf_all_samples[,, sample_iter, 1] <- sim_res[[1]]
    final_girf_all_samples[,, sample_iter, 2] <- sim_res[[2]]
    
    # for (reg_iter in 1:n_regimes) {
    #   curr_regime <- histories[[reg_iter]]
    #   regime_struc_errors <- struc_errors[[reg_iter]]
    #   nobs_reg <- nrow(curr_regime)
    #   girf_curr_reg <- array(NA, dim = c(horizon, N, n_histories))
    #   for (hist_iter in 1:n_histories) {
    #     curr_hist_pos <- sample.int(nobs_reg, size = 1)
    #     curr_hist_pos_in_ds <- pos_in_ds[[reg_iter]][curr_hist_pos]
    #     curr_hist <- curr_regime[curr_hist_pos, ] #have to take dataset[curr_hist_pos_in_ds + P, ]
    #     girf_boot <- array(NA, dim = c(horizon, N, n_boot))
    #     for (boot_iter in 1:n_boot) {
    #       curr_error_pos <- sample.int(nobs_reg, size = horizon, replace = T)
    #       curr_str_errors <- regime_struc_errors[, curr_error_pos]
    #       girf_boot[,, boot_iter] <- simul_diff_cpp(curr_B_array, N, curr_hist, curr_str_errors, struc_Sigma, 
    #                                                 horizon, reg_iter, tardata, curr_hist_pos_in_ds, 
    #                                                 P, max_delay_needed, curr_tar, tar_variable,
    #                                                 shock_var, shock_size, tar_transform)
    #     }
    #     girf_curr_reg[,, hist_iter] <- apply_girf_cpp(girf_boot)
    #   }
    #   final_girf_all_samples[,, sample_iter, reg_iter] <- apply_girf_cpp(girf_curr_reg)
    # }
    if (!quiet)
      print(paste0("Sample ", sample_iter, " of ", n_samples, " completed."), quote = F)
  }
  class(final_girf_all_samples) <- c(class(final_girf_all_samples), "GIRF_tvar")
  return(final_girf_all_samples)
}
