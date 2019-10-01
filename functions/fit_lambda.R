fit_lambda <- function(data, variables = 1, deltas = NULL, lambda_seq = seq(.1, 1, .1), P = 1, 
                       model = c("bvar", "tar", "svol"), scale = F, eps = 1e-4, ar_coeff = F, 
                       sum_of_coeff = T, tau = 10 * lambda_seq, rmse = T, model_input = NULL, plot = F) {
  # ar_p = T, stdev = T,   #taken out from function call
  if (length(tau) == 1)
    tau <- rep(tau, length(lambda_seq))
  if (length(tau) != length(lambda_seq))
    stop("Length of argument 'tau' has to be equal to length of argument 'lambda_seq' or 1.")
  model <- match.arg(model)
  data <- as.matrix(data)
  N <- ncol(data)
  N_ols <- length(variables) #variable for fit comparison
  if (is.null(deltas)) 
    deltas <- rep(0, N)
  if (scale)
    data <- scale(data)
  # if (model == "tar") {
  #   pos_tar <- which(colnames(data) == model_input$tar_variable)
  #   data_tar <- data[, variables]
  #   ols_tar <- tsDyn::TVAR(data_tar, lag = P, include = "const", model = "TAR", 
  #                          commonInter = FALSE, nthresh = 1, thDelay = 1, 
  #                          thVar = data[, pos_tar], plot = F, trace = F)
  #   rmse_ols <- colMeans(ols_tar$residuals^2)
  # } else {
  #   #OLS estimation with reduced data
  #   ols_data <- embed(data[, variables], P + 1)
  #   Y_ols <- ols_data[, 1:N_ols]
  #   X_ols <- cbind(ols_data[, -c(1:N_ols)], 1)
  #   ols_res <- .lm.fit(X_ols, Y_ols)
  #   rmse_ols <- colMeans(ols_res$residuals^2)
  # }
  ols_data <- embed(data[, variables], P + 1) #OLS on reduced dataset as benchmark
  Y_ols <- ols_data[, 1:N_ols]
  X_ols <- cbind(ols_data[, -c(1:N_ols)], 1)
  ols_res <- .lm.fit(X_ols, Y_ols)
  rmse_ols <- colMeans(ols_res$residuals^2)
  if (rmse) {
    rmse_0 <- diag(cov(diff(data)))
    rmse_ols <- rmse_ols / rmse_0[variables]
  }
  #rmse_ols <- sqrt(rmse_ols)
  #Estimating other model
  fit <- matrix(NA, length(lambda_seq), N)
  colnames(fit) <- colnames(data)
  if (model == "bvar") {
    for (i in 1:length(lambda_seq)) {
      est_res <- bvar_lambda(data = data, P = P, deltas = deltas, lambda = lambda_seq[i], tau = tau[i], 
                  eps = eps, sum_of_coeff = sum_of_coeff, train_sample = NULL)
      rmse_model <- colMeans(est_res$resids^2) #save all variables not only benchmark vars...
      if (rmse) 
        rmse_model <- rmse_model / rmse_0 #... therefore divide by rmse_0 instead of rmse_0[variables]
      #fit[i, ] <- sqrt(rmse_model)
      fit[i, ] <- rmse_model
    }
  }
  if (model == "tar") {
    for (i in 1:length(lambda_seq)) {
      est_res <- tryCatch(
        {tar_res <- bayesian_tvar(data = data, P = P, tar_variable = model_input$tar_variable, 
                                 max_d = model_input$max_d, tar_prior = model_input$tar_prior, 
                                 tar_scale = model_input$tar_scale, tar_transform = model_input$tar_transform, 
                                 reps = model_input$reps, burn = model_input$burn, 
                                 stability = model_input$stability, max_attempts = model_input$max_attempts, 
                                 deltas = deltas, lambda = lambda_seq[i], tau = tau[i], eps = eps, 
                                 sum_of_coeff = sum_of_coeff, train_sample = NULL, quiet = F, 
                                 forecast_horizon = NULL, irf_settings = NULL, 
                                 dram_settings = model_input$dram_setting)
        resids_avg <- apply(tar_res$out_resid, c(2, 3), mean)
        rmse_model <- colMeans(resids_avg^2)
        if (rmse)
          rmse_model <- rmse_model / rmse_0
        #fit[i, ] <- sqrt(rmse_model)
        fit[i, ] <- rmse_model}, 
        error = function(e) base::message("Lambda = ", lambda_seq[i], " resulted in an error."))
    }
  }
  if (model == "svol") {
    for (i in 1:length(lambda_seq)) {
      est_res <- tryCatch(
        {svol_res <- bayesian_var_svol(data = data, P = P, sigbar = model_input$sigbar, g0 = model_input$g0, 
                                      Tg0 = model_input$Tg0, reps = model_input$reps, burn = model_input$burn, 
                                      stability = model_input$stability, max_attempts = model_input$max_attempts, 
                                      deltas = deltas, lambda = lambda_seq[i], tau = tau[i], eps = eps, 
                                      sum_of_coeff = sum_of_coeff, train_sample = model_input$train_sample, 
                                      quiet = F, forecast_horizon = NULL, irf_settings = NULL)
        resids_avg <- apply(svol_res$out_resid, c(2, 3), mean)
        rmse_model <- colMeans(resids_avg^2)
        if (rmse) 
          rmse_model <- rmse_model / rmse_0
        #fit[i, ] <- sqrt(rmse_model)
        fit[i, ] <- rmse_model},
        error = function(e) base::message("Lambda = ", lambda_seq[i], " resulted in an error."))
    }
  }
  fit_variables <- apply(fit[, variables], 1, mean) #select the variables we need and then mean over them (margin = 1)
  mean_rmse_ols <-  mean(rmse_ols) #average over the variables in (small) VAR
  dist_lambda <- abs(fit_variables - mean_rmse_ols)
  min_pos <- which.min(dist_lambda)
  lambda <- lambda_seq[min_pos]
  if (plot) {
    plot(lambda_seq, dist_lambda, type = "l", ylab = "Distance", xlab = "Lambda")
    points(lambda_seq[min_pos], dist_lambda[min_pos], col = "red")
    legend("topright", paste0("Lambda: ", lambda, " , distance = ", 
                              round(dist_lambda[min_pos], 4), "."), cex = 0.6, bty = "n")
  }
  ret_list <- list(lambda = lambda, dist_lambda = dist_lambda, position = min_pos,
                   lambda_seq = lambda_seq, rmse_ols = mean_rmse_ols)
  return(ret_list)
}
