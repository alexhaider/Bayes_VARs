#only lm fit for finding lambda, no Gibbs sampling
bvar_lambda <- function(data, P, deltas = NULL, lambda = 0.1, tau = 10 * lambda, eps = 1e-4,  
                        sum_of_coeff = T, train_sample = NULL) {
  # ar_coeff = F, ar_p = T, stdev = T,  #taken out from function calls
  out_call <- match.call()
  # if (scale_vars) 
  #   data <- scale(data)
  N <- ncol(data)
  T_complete <- nrow(data) 
  # stdev <- T #for dummy prior
  # ar_p <- T
  if (is.null(train_sample))
    train_sample <- T_complete
  if (is.null(deltas))
    deltas <- rep(0, N)
  #mus <- colMeans(data[-c(1:P), ])
  dummy_p <- create_dummy_prior_man(data[1:train_sample, ], P = P, deltas = deltas, lambda = lambda, 
                                    eps = eps, sum_of_coeff = sum_of_coeff, tau = tau)
  YD <- dummy_p$Yd
  XD <- dummy_p$Xd
  data_embed <- embed(as.matrix(data), P + 1)
  T_embed <- nrow(data_embed) 
  YY <- data_embed[, 1:N] #new because remove training sample has been commented out
  XX <- cbind(data_embed[, -c(1:N)], 1) #new because remove training sample has been commented out
  Y_star <- rbind(YY, YD)
  X_star <- rbind(XX, XD)
  lm_fit <- .lm.fit(X_star, Y_star)
  resids <- lm_fit$residuals[1:T_embed, ] #no residuals from dummy obs
  ret_list <- list(out_call = out_call, resids = resids)
  return(ret_list)
}
