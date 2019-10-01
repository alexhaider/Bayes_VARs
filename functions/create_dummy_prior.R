#ar_coeff...are we using ols estimates for gamma?; otherwise 0 & 1 based on PP Test
#P...lags
#works with standard deviation or variance (stdev)
#ar_p: AR(P) or AR(1) for computing 
create_dummy_prior <- function(data, P, lambda = 0.1, tau = 10 * lambda, eps = 1e-4, ar_coeff = F, 
                               ar_p = T, stdev = T, sum_of_coeff = T) {
  N <- ncol(data)
  mus <- colMeans(data, na.rm = T)
  if (ar_p) { #same order in AR as in VAR...
    data <- embed(as.matrix(data), P + 1)
    p_ar <- P
  } else {#...or just AR(1)?
    data <- embed(as.matrix(data), 2)
    p_ar <- 1
  }
  deltas <- sigs <- pp_test <- numeric()
  for (i in 1:N) {
    lm_fit <- lm(data[, i] ~ data[, N * (1:p_ar) + i])
    deltas <- c(deltas, ifelse(coef(lm_fit)[2] > 1, 1, coef(lm_fit)[2]))
    #sigs <- c(sigs, summary(lm_fit)$sigma)
    curr_sig <- crossprod(lm_fit$residuals) / (nrow(data) - length(lm_fit$coefficients))
    sigs <- c(sigs, ifelse(stdev, sqrt(curr_sig), curr_sig))
    pp_test <- c(pp_test, suppressWarnings(tseries::pp.test(data[, i])$p.value))
  }
  pp_test_res <- pp_test
  if (!ar_coeff) {
    pp_test[pp_test > 0.05] <- 1
    pp_test[pp_test <= 0.05] <- 0
    deltas <- pp_test
  }
  Jp <- diag(1:P)
  Yd <- rbind(diag(deltas * sigs) / lambda, matrix(0, N * (P - 1), N), 
              diag(sigs), matrix(0, 1, N))
  Xd <- rbind(cbind(kronecker(Jp, diag(sigs)) / lambda, matrix(0, N * P, 1)),
              cbind(matrix(0, N, N * P), matrix(0, N, 1)),
              cbind(matrix(0, 1, N * P), eps))
  
  if (sum_of_coeff) {
    Yd <- rbind(Yd, diag(deltas * mus) / tau)
    Xd <- rbind(Xd, cbind(kronecker(matrix(1, 1, P), diag(deltas * mus) / tau), 
                          matrix(0, N, 1)))
  }
  return(list(Yd = Yd, Xd = Xd, mu = mus, delta = deltas, sig = sigs, pp_test_res = pp_test_res))
}


create_dummy_prior_man <- function(data, P, deltas = 0, lambda = 0.1, eps = 1e-4, sum_of_coeff = T, 
                                   tau = 10 * lambda) {
  N <- ncol(data)
  mus <- colMeans(data, na.rm = T)
  data <- embed(as.matrix(data), P + 1)
  if (length(deltas) == 1)
    deltas <- rep(deltas, N)
  sigs <- numeric()
  for (i in 1:N) {
    lm_fit <- lm(data[, i] ~ data[, N * (1:P) + i])
    # sigs <- c(sigs, sqrt(crossprod(lm_fit$residuals) / lm_fit$df.residual))
    sigs <- c(sigs, crossprod(lm_fit$residuals) / lm_fit$df.residual)
  }
  Jp <- diag(1:P)
  Yd <- rbind(diag(deltas * sigs) / lambda, 
              matrix(0, N * (P - 1), N), 
              diag(sigs), matrix(0, 1, N))
  Xd <- rbind(cbind(kronecker(Jp, diag(sigs)) / lambda, matrix(0, N * P, 1)),
              cbind(matrix(0, N, N * P), matrix(0, N, 1)),
              cbind(matrix(0, 1, N * P), eps))
  if (sum_of_coeff) {
    Yd <- rbind(Yd, diag(deltas * mus) / tau)
    Xd <- rbind(Xd, cbind(kronecker(matrix(1, 1, P), diag(deltas * mus) / tau), 
                          matrix(0, N, 1)))
  }
  return(list(Yd = Yd, Xd = Xd, sig = sigs))
}

# create_dummy_prior_man2 <- function(data, P, deltas = 0, lambda = 0.1, eps = 1e-4, sum_of_coeff = T,
#                                     tau = 10 * lambda, mus) {
#   N <- ncol(data)
#   #mus <- colMeans(data, na.rm = T)
#   data <- embed(as.matrix(data), P + 1)
#   if (length(deltas) == 1)
#     deltas <- rep(deltas, N)
#   sigs <- numeric()
#   for (i in 1:N) {
#     lm_fit <- lm(data[, i] ~ data[, N * (1:P) + i])
#     sigs <- c(sigs, sqrt(crossprod(lm_fit$residuals) / lm_fit$df.residual))
#   }
#   Jp <- diag(1:P)
#   Yd <- rbind(diag(deltas * sigs) / lambda, 
#               matrix(0, N * (P - 1), N), 
#               diag(sigs), matrix(0, 1, N))
#   Xd <- rbind(cbind(kronecker(Jp, diag(sigs)) / lambda, matrix(0, N * P, 1)),
#               cbind(matrix(0, N, N * P), matrix(0, N, 1)),
#               cbind(matrix(0, 1, N * P), eps))
#   if (sum_of_coeff) {
#     Yd <- rbind(Yd, diag(deltas * mus) / tau)
#     Xd <- rbind(Xd, cbind(kronecker(matrix(1, 1, P), diag(deltas * mus) / tau), 
#                           matrix(0, N, 1)))
#   }
#   return(list(Yd = Yd, Xd = Xd, sig = sigs))
# }
