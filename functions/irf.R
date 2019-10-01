# # est_res <- tvar_res1; shock_var <- 1; shock_size <- 1; type = "Triangular"
# # horizon <- 20; n_histories <- 3; n_boot <- 4; quiet <- T
# GIRF <- function(est_res, shock_var = 1, shock_size = 1, type = c("Chol", "Triangular", "Sign_restr"),
#                  n_histories = 500, n_boot = 100, horizon = 20, quiet = FALSE) {
#   #variables
#   type <- match.arg(type)
#   dataset <- eval(est_res$out_call$data) #complete dataset
#   tar_variable <- eval(est_res$out_call$tar_variable)
#   tardata <- dataset[, tar_variable]
#   P <- eval(est_res$out_call$P)
#   if (is.null(est_res$out_call$tar_transform)) {
#     tar_transform <- list(fun = "ident", inp = 1)
#   } else {
#     tar_transform <- eval(est_res$out_call$tar_transform)
#   }
#   N <- ncol(dataset)
#   YY <- est_res$model_specific$data_embed[, 1:N]
#   XX <- cbind(est_res$model_specific$data_embed[, -c(1:N)], 1)
#   ZZ <- est_res$model_specific$ZZ
#   horizon <- horizon + 1
#   n_regimes <- 2
#   n_samples <- est_res$out_call$reps - est_res$out_call$burn
#   
#   if (is.character(shock_var)) {
#     char_shock_var <- shock_var
#     shock_var <- which(colnames(dataset) == shock_var)
#     if (length(shock_var) == 0)
#       stop("Error. Variable ", char_shock_var, " does not exist!")
#   }
#   if (!quiet) {
#     cat(("GIRF\n====\nshocked variable: "), colnames(dataset)[shock_var], 
#         "\nshock size:", shock_size, "\nGirf Simulation:\n")
#     count_prog <- 0
#   }
#   
#   #Loop over samples
#   sample_iter <- 1
#   final_girf_all_samples <- array(NA, dim = c(horizon, N, n_regimes, n_samples))
#   for (sample_iter in 1:n_samples) {
#     curr_B1 <- create_Bmat(est_res$out_beta1[sample_iter, ], N)
#     curr_B2 <- create_Bmat(est_res$out_beta2[sample_iter, ], N)
#     curr_B_array <- array(c(curr_B1, curr_B2), dim = c(nrow(curr_B1), ncol(curr_B1), n_regimes))
#     curr_delay <- est_res$out_delay[sample_iter]
#     max_delay_needed <- curr_delay + tar_transform$inp - 1
#     curr_tar <- est_res$out_tar[sample_iter]
#     reg1 <- ZZ[, curr_delay] <= curr_tar; reg2 <- !reg1
#     pos_in_ds <- list(reg1 = which(reg1), reg2 = which(reg2)) #saving postiiotns of obs in ds
#     XX_1 <- XX[reg1,, drop = F]; XX_2 <- XX[reg2,, drop = F]
#     resid_reg1 <- YY[reg1,, drop = F] - XX_1 %*% t(curr_B1)
#     resid_reg2 <- YY[reg2,, drop = F] - XX_2 %*% t(curr_B2)
# 
#     #Prepare Sigma, histories and structural errors
#     struc_Sigma <- struc_errors <- histories <- list()
#     struc_Sigma[[1]] <- struc_covmat_cpp(est_res$out_sigma1[sample_iter,, ], type) #lower triang
#     struc_Sigma[[2]] <- struc_covmat_cpp(est_res$out_sigma2[sample_iter,, ], type)
#     struc_errors[[1]] <- struc_error_cpp(struc_Sigma[[1]], resid_reg1)
#     struc_errors[[2]] <- struc_error_cpp(struc_Sigma[[2]], resid_reg2)
#     histories[[1]] <- XX_1; histories[[2]] <- XX_2
#     
#     # reg_iter <- hist_iter <- boot_iter <- 1
#     final_girf <- array(NA, dim = c(horizon, N, n_regimes))
#     #girf_by_reg <- list()
#     for (reg_iter in 1:n_regimes) {
#       curr_regime <- histories[[reg_iter]]
#       regime_struc_errors <- struc_errors[[reg_iter]]
#       hist_pos <- sample(1:nrow(curr_regime), size = n_histories, replace = T)
#       error_pos <- sample(1:ncol(regime_struc_errors), size = n_histories * n_boot * horizon, replace = T)
#       girf_curr_reg <- array(NA, dim = c(horizon, N, n_histories))
#       for (hist_iter in 1:n_histories) {
#         curr_hist_pos <- hist_pos[hist_iter]
#         curr_hist_pos_in_ds <- pos_in_ds[[reg_iter]][curr_hist_pos]
#         curr_hist <- curr_regime[curr_hist_pos, ]
#         girf_boot <- array(NA, dim = c(horizon, N, n_boot))
#         for (boot_iter in 1:n_boot) {
#           curr_error_pos <- (hist_iter - 1) * (boot_iter * horizon - 1) + 
#             (boot_iter - 1) * horizon
#           curr_error_pos <- (curr_error_pos + 1):(curr_error_pos + horizon) #position of errors
#           curr_str_errors <- regime_struc_errors[, error_pos[curr_error_pos], drop = F] #get errors
#           curr_str_errors_add_shock <- curr_str_errors
#           curr_str_errors_add_shock[shock_var, 1] <- shock_size
#           reduced_form_errors <- struc_Sigma[[reg_iter]] %*% curr_str_errors
#           reduced_form_errors_add_shock <- struc_Sigma[[reg_iter]] %*% curr_str_errors_add_shock
#           tard_short <- tardata[(curr_hist_pos_in_ds + P - max_delay_needed):(curr_hist_pos_in_ds + P - 1)]
#           # simul_wo_shock <- simul_sys(curr_B_array, N, curr_hist, reduced_form_errors, horizon, reg_iter,
#           #                             tard_short, curr_tar, tar_transform, tar_variable)
#           # sim_res <- simul_diff_cpp(curr_B_array, N, curr_hist, reduced_form_errors, reduced_form_errors_add_shock,
#           #                          horizon, reg_iter, tard_short, curr_tar, tar_transform, tar_variable)
#           # simul_wo_shock <- simul_sys_cpp(curr_B_array, N, curr_hist, reduced_form_errors, horizon, reg_iter,
#           #                                 tard_short, curr_tar, tar_transform, tar_variable)
#           # simul_add_shock <- simul_sys_cpp(curr_B_array, N, curr_hist, reduced_form_errors_add_shock, 
#           #                                  horizon, reg_iter, tard_short, curr_tar, 
#           #                                  tar_transform, tar_variable)
#           #  simul_add_shock <- simul_sys(curr_B_array, N, curr_hist, reduced_form_errors_add_shock, 
#           #                               horizon, reg_iter, tard_short, curr_tar, 
#           #                               tar_transform, tar_variable)
#           girf_boot[,, boot_iter] <- simul_diff_cpp(curr_B_array, N, curr_hist, reduced_form_errors, 
#                                                     reduced_form_errors_add_shock,
#                                                     horizon, reg_iter, tard_short, 
#                                                     curr_tar, tar_transform, tar_variable)
#           if (!quiet) {
#             count_prog <- count_prog + 1
#             prog <- count_prog / (n_boot * n_histories * n_regimes * n_samples) * 100
#             if (prog %% 1 == 0)
#               cat(prog)
#           }
#         }
#         girf_curr_reg[,, hist_iter] <- apply_girf_cpp(girf_boot)
#         #girf_curr_reg[,, hist_iter] <- apply(girf_boot, MARGIN = c(1, 2), mean)
#         #girf_by_reg[[reg_iter]] <- girf_curr_reg
#       }
#       final_girf[,, reg_iter] <- apply_girf_cpp(girf_curr_reg)
#       #final_girf[,, reg_iter] <- apply(girf_curr_reg, MARGIN = c(1, 2), mean)
#     }
#     final_girf_all_samples[,,, sample_iter] <- final_girf
#   }
#   #ret_list <- list(girf = final_girf, girf_by_reg = girf_by_reg, struc_Sigma = struc_Sigma)
#   class(final_girf_all_samples) <- c(class(final_girf_all_samples), "GIRF_tvar")
#   return(final_girf_all_samples)
# }
# 
# #history <- curr_hist; shock_sequence <- reduced_form_errors; start_reg <- reg_iter; tar_transf = tar_transform
# simul_sys <- function(curr_B_array, N, history, shock_sequence, horizon, start_reg, 
#                       tard_short, curr_tar, tar_transf, tar_variable) {
#   #history + 1 because we evaluate next regime at the end of the loop which means it is 
#   #evaluated once for nothing
#   curr_B <- curr_B_array[,, start_reg]
#   len_history <- length(history) - 1 #makes it easier to delete old history values later on
#   sim_res <- matrix(NA, nrow = horizon, ncol = N)
#   for (i in 1:horizon) {
#     sim_res[i, ] <- curr_B %*% history + shock_sequence[, i, drop = F]
#     #history <- c(sim_res[i, ], history)
#     tard_short <- c(tard_short, sim_res[i, tar_variable])
#     ZZ_used <- transform_trans_variable(tard_short, tar_transf)
#     new_reg <- as.numeric(ZZ_used[i + 1] > curr_tar) + 1 #which regime is next
#     #ind_fn_history[i + 1, new_reg] <- 1
#     curr_B <- curr_B_array[,, new_reg] #update coefficient matrix
#     history <- c(c(sim_res[i, ], history)[1:len_history], 1)
#   }
#   return(sim_res)
# }
# 
# transform_trans_variable <- function(x, tar_transform) {
#   #tar_transform <- eval(tar_transform)
#   if (tar_transform$fun == "ident")
#     return(x)
#   if (tar_transform$fun == "growth") 
#     return(log(x) - log(dplyr::lag(x, tar_transform$inp)))
#   if (tar_transform$fun == "MA")
#     return(zoo::rollmean(x, tar_transform$inp, fill = NA, align = "right"))
#   if (tar_transform$fun == "sum")
#     return(zoo::rollsum(x, tar_transform$inp, fill = NA, align = "right"))
# }

irf <- function(type = c("Chol", "Triangular", "Sign_restr"), B_mat, Sigma, shocked_variable, 
                shock_size, horizon, restrict = NULL) {
  N <- nrow(B_mat)
  P <- (ncol(B_mat) - 1) / N
  type <- match.arg(type)
  if (type == "Chol") 
    A0 <- chol(Sigma) #upper triang
  if (type == "Triangular") { #does it really make sense?
    chol_sigma <- t(chol(Sigma))
    Gamma_power_minus_onehalf <- diag(1 / diag(chol_sigma))
    A0 <- t(chol_sigma %*% Gamma_power_minus_onehalf) #has to be transposed because row 17
  }
  yhat <- matrix(0, nrow = P + horizon, ncol = N)
  v <- rep(0, N)
  v[shocked_variable] <- shock_size
  yhat_index <- P:1
  yhat[P + 1, ] <- c(t(yhat[yhat_index, ]), 0) %*% t(B_mat) + v %*% A0 #beacuse transp
  yhat_index <- yhat_index + 1
  #intercept = 0 -> doesn't matter in IR
  for (j in (P + 2):nrow(yhat)) {
    yhat[j, ]  <- c(t(yhat[yhat_index, ]), 0) %*% t(B_mat)
    yhat_index <- yhat_index + 1
  }
  return(yhat[-c(1:P), ])
}