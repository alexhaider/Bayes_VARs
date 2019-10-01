AR_sv <- function(Y, p = 4, intercept = TRUE,  h = 4, ar_priors, sv_priors, b_start = NULL, 
                  h_start = NULL, reps = 10000, burn = 5000) { #works compared with stochvol package fc
  #extract prior
  b0 <- ar_priors$b0
  B0inv <- solve(ar_priors$B0)
  priormu <- sv_priors$priormu
  priorphi <- sv_priors$priorphi
  priorsigma <- sv_priors$priorsigma
  #prepare data
  y_emb <- embed(Y, p + 1)
  y <- y_emb[, 1]
  if (intercept) {
    X <- cbind(1, y_emb[, -1]) 
    colnames(X) <- c("const", paste0("lag", 1:p))
  } else {
    X <- y_emb[, -1]
    colnames(X) <- paste0("lag", 1:p)
  } 
  N <- length(y)
  #starting vals
  len_b <- length(b0)
  if (is.null(b_start))
    b_start <- rep(0, len_b) #b0 is the prior
  b_vec <- b_start
  if (is.null(h_start)) 
    h_start <- list(para = c(mu = 0, phi = .9, sigma = .1), latent = rep(0, N)) #starting
  h_new <- h_start
  #saving results
  draws <- matrix(NA_real_, nrow = reps - burn, ncol = 3 + len_b)
  colnames(draws) <- c("mu", "phi", "sigma", paste("beta", 0:(len_b - 1), sep = "_"))
  draws_latent <- matrix(NA_real_, nrow = N, ncol = reps - burn)
  colnames(draws_latent) <- paste("fc", 1:(reps - burn), sep = "_")
  y_predict <- matrix(NA_real_, nrow = h + p, ncol = reps - burn)
  y_predict[1:p, ] <- tail(y, p)
  colnames(y_predict) <- paste("fc", 1:(reps - burn), sep = "_")
  for (i in 1:reps) {
    eps <- y - X %*% b_vec #get residuals for stochvol sampler
    eps[which(eps == 0)] <- 1e-4 #no zeros for h_new
    h_new <- svsample2(eps, startpara = para(h_new), startlatent = latent(h_new), 
                       priormu = priormu, priorphi = priorphi, priorsigma = priorsigma)
    normalizer <- as.numeric(exp(-latent(h_new) / 2)) #heterosked remover
    X_new <- X * normalizer #remove heterosked
    y_new <- y * normalizer
    # Sigma <- solve(crossprod(X_new) + B0inv)
    Sigma <- xpx_inv_offset_cpp(X_new, B0inv) #(X'X + B0inv)^-1
    mu_ar <- Sigma %*% (crossprod(X_new, y_new) + B0inv %*% b0)
    b_vec <- as.numeric(mvrnorm_cpp(1, mu_ar, Sigma))
    #saving results after burn in
    if (i > burn) {
      params_sv <- para(h_new) #latest sampled parameters for state equation
      latent_sv <- latent(h_new) #sampled h_s
      draws[i - burn, 1:3]  <- params_sv
      draws[i - burn, 4:(len_b + 3)] <- b_vec
      draws_latent[, i - burn] <- latent_sv
      # draws[i - burn, (len_b + 4):ncol(draws)] <- latent(h_new)
      #forecast
      h_predict <- rep(NA_real_, h + 1)
      h_predict[1] <- latent_sv[N]
      for (j in 1:h) {
        h_predict[j + 1] <- params_sv[1] + params_sv[2] * (h_predict[j] - params_sv[1]) +
          rnorm(1, 0, params_sv[3])
        if (intercept) {
          y_predict[j + p, i - burn] <- sum(c(1, rev(y_predict[j:(j + p - 1) , i - burn])) * b_vec) +
            rnorm(1, 0, exp(h_predict[j + 1] / 2))
        } else {
          y_predict[j + p, i - burn] <- sum(rev(y_predict[j:(j + p - 1) , i - burn]) * b_vec) +
            rnorm(1, 0, exp(h_predict[j + 1] / 2))
        }
      }
    }
    if (i %% 5000 == 0)
      print(paste0("SVOL_AR, Replication ", i, " of ", reps, "."), quote = F)
  }
  ret_list <- list(draws = draws, draws_latent = draws_latent, y_predict = y_predict,
                   dataset = Y)
  return(ret_list)
}

AR_sv_rw <- function(Y, p = 4, intercept = TRUE,  h = 4, ar_priors, sv_priors, b_start = NULL, 
                     h_start = NULL, reps = 10000, burn = 5000) {
  #extract prior
  b0 <- ar_priors$b0
  B0inv <- solve(ar_priors$B0)
  g0 <- sv_priors$g0
  Tg0 <- sv_priors$Tg0
  mubar <- sv_priors$mubar
  sigbar <- sv_priors$sigbar
  #prepare data
  y_emb <- embed(Y, p + 1)
  y <- y_emb[, 1]
  if (intercept) {
    X <- cbind(1, y_emb[, -1]) 
    colnames(X) <- c("const", paste0("lag", 1:p))
  } else {
    X <- y_emb[, -1]
    colnames(X) <- paste0("lag", 1:p)
  } 
  N <- length(y)
  len_b <- length(b0)
  #starting vals
  if (is.null(b_start))
    b_start <- rep(0, len_b) #b0 is the prior
  b_vec <- b_start
  if (is.null(h_start))  {
    h_start <- diff(y)^2
    h_start <- c(h_start[1], h_start, h_start[length(h_start)]) + sd(y) / 1000
  }
  # h_start <- rep(1, N + 1) #starting
  h_new <- as.matrix(h_start, ncol = 1)
  #starting loop
  draws <- matrix(NA_real_, nrow = reps - burn, ncol = 1 + len_b) #"1 +" for g
  colnames(draws) <- c("g", paste("beta", 0:(len_b - 1), sep = "_"))
  draws_latent <- matrix(NA_real_, nrow = N, ncol = reps - burn)
  colnames(draws_latent) <- paste("fc", 1:(reps - burn), sep = "_")
  y_predict <- matrix(NA_real_, nrow = h + p, ncol = reps - burn)
  y_predict[1:p, ] <- tail(y, p) 
  colnames(y_predict) <- paste("fc", 1:(reps - burn), sep = "_")
  for (i in 1:reps) {
    eps <- y - X %*% b_vec
    # eps[which(eps == 0)] <- 1e-4 #no zeros for h_new
    g <- sample_g(h_new, Tg0, g0)
    h_new <- svol_mh(h_new, eps, g, mubar, sigbar)
    # h_new <- h_new[-1] #getting rid of starting value
    normalizer <- as.numeric(exp(-log(h_new[-1]) / 2)) #heterosked remover
    X_new <- X * normalizer
    y_new <- y * normalizer
    Sigma <- xpx_inv_offset_cpp(X_new, B0inv)
    mu_ar <- Sigma %*% (crossprod(X_new, y_new) + B0inv %*% b0)
    b_vec <- as.numeric(mvrnorm_cpp(1, mu_ar, Sigma))
    #saving results after burn in
    if (i > burn) {
      draws[i - burn, 1]  <- g
      draws[i - burn, 2:(len_b + 1)] <- b_vec
      draws_latent[, i - burn] <- h_new[-1]
      # draws[i - burn, (len_b + 4):ncol(draws)] <- latent(h_new)
      #forecast
      h_predict <- rep(NA_real_, h + 1)
      h_predict[1] <- log(h_new[N + 1])
      for (j in 1:h) {
        h_predict[j + 1] <- h_predict[j] + rnorm(1, 0, sqrt(g)) #g or sqrt(g)
        if (intercept) {
          y_predict[j + p, i - burn] <- sum(c(1, rev(y_predict[j:(j + p - 1) , i - burn])) * b_vec) +
            rnorm(1, 0, exp(h_predict[j + 1] / 2))
        } else {
          y_predict[j + p, i - burn] <- sum(rev(y_predict[j:(j + p - 1) , i - burn]) * b_vec) +
            rnorm(1, 0, exp(h_predict[j + 1] / 2))
        }
      }
    }
    if (i %% 5000 == 0)
      print(paste0("SVOL_RW, Replication ", i, " of ", reps, "."), quote = F)
  }
  ret_list <- list(draws = draws, draws_latent = draws_latent, y_predict = y_predict,
                   dataset = Y)
  return(ret_list)
}

AR_hs <- function(Y, p = 4, intercept = TRUE, h = 4, ar_priors, sigma_start = NULL,
                  h_start = NULL, reps = 10000, burn = 5000) {
  #extract prior
  b0 <- ar_priors$b0
  B0inv <- solve(ar_priors$B0)
  c0 <- ar_priors$c0
  C0 <- ar_priors$C0
  #prepare data
  y_emb <- embed(Y, p + 1)
  y <- y_emb[, 1]
  if (intercept) {
    X <- cbind(1, y_emb[, -1]) 
    colnames(X) <- c("const", paste0("lag", 1:p))
  } else {
    X <- y_emb[, -1]
    colnames(X) <- paste0("lag", 1:p)
  } 
  N <- length(y)
  len_b <- length(b0)
  if (is.null(sigma_start))
    sigma_start <- 1
  sigma2_draw <- sigma_start #starting value for sigma2, no starting value for beta needed
  #saving res
  draws <- matrix(NA_real_, nrow = reps - burn, ncol = 1 + len_b) #"1 +" for sigma
  colnames(draws) <- c( "sigma", paste("beta", 0:(len_b - 1), sep = "_"))
  y_predict <- matrix(NA_real_, nrow = h + p, ncol = reps - burn)
  colnames(y_predict) <- paste("fc", 1:(reps - burn), sep = "_")
  y_predict[1:p, ] <- tail(y, p) 
  #pre-calcs 
  preCov <- solve(crossprod(X) + B0inv)
  preMean <- preCov %*% (crossprod(X, y) + B0inv %*% b0)
  preDf <- c0 + N / 2 + p / 2
  for (i in 1:reps) {
    b_vec <- as.numeric(mvrnorm_cpp(1, preMean, sigma2_draw * preCov))
    tmp <- C0 + .5 * (crossprod(y - X %*% b_vec) + 
                        crossprod((b_vec - b0), B0inv) %*% (b_vec - b0))
    sigma2_draw <- 1 / rgamma(1, preDf, rate = tmp) #use rate!!!
    #saving results after burn in
    if (i > burn) {
      draws[i - burn, ]  <- c(sqrt(sigma2_draw), b_vec)
      for (j in 1:h) {
        if (intercept) {
          y_predict[j + p, i - burn] <- sum(c(1, rev(y_predict[j:(j + p - 1) , i - burn])) * b_vec) +
            rnorm(1, 0, sqrt(sigma2_draw))
        } else {
          y_predict[j + p, i - burn] <- sum(rev(y_predict[j:(j + p - 1) , i - burn]) * b_vec) +
            rnorm(1, 0, sqrt(sigma2_draw))
        }
      }
    }
    if (i %% 5000 == 0)
      print(paste0("AR, Replication ", i, " of ", reps, "."), quote = F)
  }
  ret_list <- list(draws = draws, y_predict = y_predict, dataset = Y)
  return(ret_list)
}

univ_star <- function(Y, p = 4, intercept = TRUE, h = 4, ar_priors, st_priors, d0 = 1,
                    delta_gamma = 0.02, delta_cc = 0.02, #for drawing new params MH
                    sigma_start = NULL, reps = 10000, burn = 5000) {
  zero_vec <- as.matrix(0) #so I can use lik_cpp which is multivariate and needs a matrix
  b0 <- ar_priors$b0 #prior for beta
  B0inv <- solve(ar_priors$B0)
  c0 <- ar_priors$c0 #prior for sigma
  C0 <- ar_priors$C0
  smooth_a <- st_priors$smooth_a #prior for gamma (smoothness)
  smooth_b <- st_priors$smooth_b
  thresh_m <- st_priors$thresh_m #prior for center 
  thresh_sd <- sqrt(st_priors$thresh_v) #To standard dev
  #prepare data
  y_emb <- embed(Y, p + 1)
  yy <- y_emb[, 1]
  if (intercept) {
    xx <- cbind(1, y_emb[, -1]) 
    colnames(xx) <- c("const", paste0("lag", 1:p))
  } else {
    xx <- y_emb[, -1]
    colnames(xx) <- paste0("lag", 1:p)
  } 
  N <- length(yy)
  len_b <- length(b0)
  #starting vals
  if (is.null(sigma_start))
    sigma_start <- 1
  sigma2_draw <- sigma_start
  gamma <- rgamma(1, smooth_a, smooth_b)
  # cc <- rnorm(1, thresh_m, thresh_sd)
  cc <- mean(yy) #starting at mean; simple
  d <- sample(1:d0, 1)
  #saving res
  draws <- matrix(NA_real_, nrow = reps - burn, ncol = 4 + len_b) 
  colnames(draws) <- c("gamma", "thresh", "delay", "sigma", 
                       paste("beta", 0:(len_b - 1), sep = "_"))
  y_predict <- matrix(NA_real_, nrow = h + p, ncol = reps - burn)
  colnames(y_predict) <- paste("fc", 1:(reps - burn), sep = "_")
  y_predict[1:p, ] <- tail(yy, p) 
  #pre-calcs 
  preDf <- c0 + N / 2 + p / 2
  n_accept <- 0
  for (i in 1:reps) {
    ss <- xx[, d + 1] #threshold observations
    G <- plogis(ss, cc, 1 / gamma) #transition function
    zz <- cbind(xx, xx * G) #regressors
    #sample b
    Sigma <- solve(B0inv + crossprod(zz))
    mu_ar <- Sigma %*% (crossprod(zz, yy) + B0inv %*% b0)
    b_vec <- as.numeric(mvrnorm_cpp(1, mu_ar, sigma2_draw * Sigma))
    #sample sigma
    resid <- yy - zz %*% b_vec
    tmp <- C0 + .5 * (crossprod(resid) + crossprod((b_vec - b0), B0inv) %*% (b_vec - b0))
    sigma2_draw <- 1 / rgamma(1, preDf, rate = tmp)
    #sample gamma, cc
    gamma_star <- rgamma(1, gamma^2 / delta_gamma, gamma / delta_gamma)
    cc_star <- rnorm(1, cc, delta_cc)
    G_new <- plogis(ss, cc_star, 1 / gamma_star)
    z_new <- cbind(xx, xx * G_new)
    resid_new <- yy - z_new %*% b_vec
    lh_new <- lik_cpp(resid_new, zero_vec, as.matrix(sigma2_draw))
    prior_new <- dgamma(gamma_star, smooth_a, smooth_b, log = T) + 
      dnorm(cc_star, thresh_m, thresh_sd, log = T)
    lh_old <- lik_cpp(resid, zero_vec, as.matrix(sigma2_draw))
    prior_old <- dgamma(gamma, smooth_a, smooth_b, log = T) + 
      dnorm(cc, thresh_m, thresh_sd, log = T)
    if (exp(lh_new + prior_new - lh_old - prior_old) > runif(1)) {
      gamma <- gamma_star
      cc <- cc_star
      n_accept <- n_accept + 1
    }
    if (i > 100 && i %% 50 == 0 && i < (burn - 100)) {
      if (n_accept / i > .5) {
        delta_gamma <- delta_gamma * 1.01
        delta_cc <- delta_cc * 1.01
      } 
      if (n_accept / i < .15) {
        delta_gamma <- delta_gamma * .99
        delta_cc <- delta_cc * .99
      }
    }
    probs_d <- rep(NA, d0)
    for (j in 1:d0) {
      ss <- xx[, j + 1]
      G <- plogis(ss, cc, 1 / gamma)
      z <- cbind(xx, xx * G)
      resid <- yy - z %*% b_vec
      probs_d[j] <- lik_cpp(resid, zero_vec, as.matrix(sigma2_draw))
    }
    probs_d <- exp(probs_d - max(probs_d))
    d <- which(rmultinom(n = 1, size = 1, prob = probs_d) == 1)
    #saving results after burn in
    if (i > burn) {
      draws[i - burn, ]  <- c(gamma, cc, d, sqrt(sigma2_draw), b_vec)
      ss_predict <- rep(NA, h) #threshold variable
      ss_predict[1:d] <- tail(yy, d) #last y value or second to last y value? should be last i think
      xx_predict <- c(tail(yy, p), rep(NA, h)) #saving the last y values and predicted ys
      for (j in 1:h) {
        G_predict <- plogis(ss_predict[j], cc, 1 / gamma)
        if (intercept) {
          zz_predict <- c(1, rev(xx_predict[j:(j + p - 1)]), G_predict, #intercept, values, intercept,..
                          rev(xx_predict[j:(j + p - 1)] * G_predict))
          tmp <- sum(zz_predict * b_vec) + rnorm(1, 0, sqrt(sigma2_draw))
          y_predict[j + p, i - burn] <- tmp
        } else {
          zz_predict <- c(rev(xx_predict[j:(j + p - 1)]),
                          rev(xx_predict[j:(j + p - 1)] * G_predict))
          tmp <- sum(zz_predict * b_vec) + rnorm(1, 0, sqrt(sigma2_draw))
          y_predict[j + p, i - burn] <- tmp
        }
        xx_predict[j + p] <- tmp
        if (j + d - 1 < length(ss_predict))
          ss_predict[j + d] <- tmp #add to delayed thresh
      }
    }
    if (i %% 5000 == 0)
      print(paste0("STAR, Replication ", i, " of ", reps, ". Accept ratio: ",  
                   round(n_accept / i, 5), "."), quote = F)
  }
  ret_list <- list(draws = draws, y_predict = y_predict, dataset = Y)
  return(ret_list)
}

