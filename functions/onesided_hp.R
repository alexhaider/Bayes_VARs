hp_onesided_kalman <- function(y, lambda = 1600, initial_val = NULL, initial_mse = NULL) {
  #Kalman filter to optimally one-sidedly filter the series that renders the standard two-sided HP filter optimal.
  #       Input: y - a Txn data matrix, where T is the number of observations on n variables
  #              lambda - a scalar. This is a smoothing parameter.
  #              initial_val - a 2xn matrix with initial values of the state estimate for each variable in y.
  #                            The underlying state vector is 2x1. Hence two values are needed for each variable in y.
  #                            Optional: if not entered, default backwards extrapolations based on the
  #                            first two observations will be used.
  #              initial_mse - a list with n elements, each element being a 2x2 matrix of initial MSE estimates
  #                            Optional: if not entered, default matrix with large variances used.
  #       Output: y_trend - a Txn matrix of extracted trends for each of the n variables.
  #               y_cycle - a Txn matrix of deviations from the extracted trends for each of the n variables.
  #
  #       Stock J.H. and M.W. Watson (1999). "Forecasting inflation," Journal of Monetary Economics,
  #                                      vol. 44(2), pages 293-335, October.
  #       The one-sided HP trend estimate is constructed as the Kalman
  #       filter estimate of tau_t in the model:
  #       y_t=tau_t+epsilon_t
  #       (1-L)^2 tau_t=eta_t"
  if (is.null(dim(y))) {
    y <- matrix(y, ncol = 1)
  } else {
    y <- matrix(y, ncol = ncol(y))
  }
  K <- ncol(y)
  nobs <- nrow(y)
  if (!is.null(initial_val)) {
    if (!is.matrix(initial_val) || dim(initial_val) != c(2, K))
      stop("'initival_val' has to be NULL or a matrix of dimension ", 2, "x", K)
  }
  if (!is.null(initial_mse)) {
    if (is.list(initial_mse)) {
      if (length(initial_mse) != K || any(unique(as.numeric(sapply(initial_mse, dim))) != 2))
        stop("'initial_mse' has to be a list of length ", K ,". All elements of 'initial_mse' have
             to be matrices of dimension 2x2.")
    } else {
      stop("'initial_mse' has to be NULL or a list.")
    }
  }
  q <- 1 / lambda #the signal-to-noise ration: i.e. var eta_t / var epsilon_t
  F_mat <- matrix(c(2, 1, -1, 0), nrow = 2) #(1 - L)^2 \tau_t = \eta_t; The state transition matrix
  H_mat <- matrix(c(1, 0), nrow = 1) #The observation matrix
  Q_mat <- matrix(c(q, 0, 0, 0), nrow = 2) #The variance-covariance matrix of the errors in the state equation
  R <- 1 #%The variance of the error in the observation equation

  y_trend <- matrix(NA, nrow = nobs, ncol = K)
  for (k in 1:K) {
    if (is.null(initial_val)) {
      # If the user didn't provide an intial value for state estimate,
      # extrapolate back two periods from the observations
      x_val <- matrix(c(2 * y[1, k] - y[2, k], 3 * y[1, k] - 2 * y[2, k]), nrow = 2)
    } else {
      x_val <- initial_val[, k]
    }
    if (is.null(initial_mse)) {
      # If the user didn't provide an intial value for the MSE, set a rather high one
      P_mat <- matrix(c(1e5, 0, 0, 1e5), nrow = 2)
    } else {
      P_mat <- initial_mse[[k]]
    }
    for (iter in 1:nobs) {
      #Kalman Update
      S_mat <- H_mat %*% P_mat %*% t(H_mat) + R
      K_mat <- F_mat %*% P_mat %*% t(H_mat) %*% solve(S_mat)
      x_val <- F_mat %*% x_val + K_mat %*% (y[iter, k] - H_mat %*% x_val) #State estimate
      temp <- F_mat - K_mat %*% H_mat
      P_mat <- temp %*% P_mat %*% t(temp)
      P_mat <- P_mat + Q_mat + K_mat %*% R %*% t(K_mat) #MSE
      #save new trend value
      y_trend[iter, k] <- x_val[2]
    }
  }
  y_cycle <- y - y_trend
  ret_list <- list(trend = y_trend, cycle = y_cycle)
  return(ret_list)
}

