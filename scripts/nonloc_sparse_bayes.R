library(tidyverse)
library(LaplacesDemon)
library(FactoMineR)
library(MatrixCorrelation)
source("MH.R")
source("pmom.R")
source("pos_eta.R")
source("eta_init.R")
get_mode <- function(x){
  # function  used for obtaining mode from a series of numbers
  uniq_vals <- unique(x)
  freq <- tabulate(match(x, uniq_vals))
  uniq_vals[which.max(freq)]
}
trim <- function(x){
  x <- as.matrix(x)
  idx.tail <- apply(x, 2, function(col){any(col!=0)}) %>%
    which %>% max
  return(as.matrix(x[,1:idx.tail]))
}

################### Sparse Bayesian Infinite Factor Model ####################
nonloc_sparse_bayes <- function(Y, eta0=NULL, 
                                num_iter=500, 
                                num_burn=200, thin=5, 
                                num_slice=(num_iter-num_burn)/thin){
  ### Parameter initialization ###
  N <- nrow(Y)
  p <- ncol(Y)
  std <- apply(Y, 2, sd) # std of Y
  # Y <- scale(Y)
  
  #k_tilde <- round(5*log(p)) # initial guess of k
  k_tilde = 3
  k_ast <- k_tilde
  
  psi <- 0.01 # dispersion parameter for pMOM density
  p0 <- 0.1
  a_p0 = b_p0 <- 5 # hyperparameter for the prior of p0
  a_sigma <- 1
  b_sigma <- 0.3
  a1 <- 2.1 # shape parameter for \delta_1 prior
  b1 <- 1 # rate parameter for \delta_1 prior
  a2 <- 3.1 # shape parameter for \delta_l prior
  b2 <- 1 # rate parameter for \delta_l prior
  df <- 3 # degree of freedom
  alpha0 <- -1 # \alpha_0 in adaptation probability
  alpha1 <- -5e-4 # \alpha_1 in adaptation probability
  epsilon <- 1e-4 # threshold for setting elements to 0
  
  Z = eta <- matrix(0, N, k_ast)
  Lambda <- matrix(0, nrow=p, ncol=k_ast)
  Sigma_inv <- rgamma(p, shape=a_sigma, rate=b_sigma)
  Sigma <- diag(1/Sigma_inv)
  phi <- rgamma(p*k_ast, shape=df/2, rate=df/2) %>% matrix(nrow=p, ncol=k_ast)
  delta <- c(rgamma(1, shape=a1, rate=b1), rgamma(k_ast-1, shape=a2, rate=b2))
  tau <- cumprod(delta)
  
  # A cube to store eta in each iteration after the burn-in phase
  #cube_eta <- array(0, dim=c(N, k_ast, num_slice))
  
  #k_est <- rep(0, num_iter) # list to record k_ast in each iteration
  
  # saving RV
  RV_Z <- rep(0, num_iter)
  RV_L <- rep(0, num_iter)
  RV_L_2 <- rep(0, num_iter)
  RV_eta <- rep(0, num_iter)
  RV_eta_2 <- rep(0, num_iter)
  RV_L_eta <- rep(0, num_iter)
  
  
  ### Gibbs sampling ###
  for (iter in 1:num_iter){
    # update Z
    Z <- sapply(1:N, function(i){
      sapply(1:k_ast, function(k){
        a <- 1 / psi + Lambda[, k] %*% diag(Sigma_inv) %*% Lambda[, k]
        b <- t(Y[i,] - Lambda[, -k] %*% eta[i, -k]) %*% diag(Sigma_inv) %*% Lambda[, k]
        t <- (a * psi)^(-3/2) * (1 + b^2 / a) * exp(b^2 / (2 * a))
        prob <- 1 - (1 - p0) / (1 - p0 + p0 * t)
        # prob
        rbinom(1, 1, prob)
      })
    }) %>% t
    
    # update eta
    idx_slab <- which(Z==1, arr.ind=T)
    idx_spike <- which(Z==0, arr.ind=T)
    
    eta_slab <- sapply(1:nrow(idx_slab), function(m){
      i <- idx_slab[m, 1]
      k <- idx_slab[m, 2]
      a <- 0.5 * (Sigma_inv %*% Lambda[, k]^2 + 1/psi)
      b <- 0.5 * t(Lambda[, k]) %*% diag(Sigma_inv) %*% (Lambda[, -k] %*% eta[i, -k] - Y[i,])
      MH(init=0.01, func=pos_eta, a=a, b=b)
    })
    eta[idx_slab] <- eta_slab
    eta[idx_spike] <- 0
    
    # update p0
    p0 <- rbeta(1, sum(sum(Z))+a_p0, N*k_ast-sum(sum(Z))+b_p0)
    
    # update Lambda
    Lambda <- sapply(1:p, function(j){
      D_inv <- diag(phi[j,] * tau)
      # covariance of lambda posterior
      cov_lambda <- (D_inv + Sigma_inv[j] * t(eta) %*% eta) %>%
        chol %>% chol2inv
      # mean of lambda posterior
      mean_lambda <- (Sigma_inv[j] * cov_lambda %*% t(eta) %*% Y[,j]) %>%
        as.vector
      rmvn(1, mean_lambda, cov_lambda)
    }) %>% t
    
    # update phi
    rate_phi <- df/2 + Lambda^2 %*% diag(tau) # rate of phi posterior
    phi <- sapply(1:p, function(j){
      sapply(1:k_ast, function(h){
        rgamma(1, shape=(df+1)/2, rate=rate_phi[j,h])
      })
    }) %>% t
    
    # update delta & tau
    phiLambda <- apply(phi*Lambda^2, 2, sum) # the product of phi' and Lambda
    for (h in 1:k_ast){
      if (h == 1){
        ad <- a1 + p * k_ast / 2
      } else {
        ad <- a2 + p * (k_ast - h + 1) / 2
      }
      bd <- 1 + 0.5 / delta[h] * tau[h:k_ast] %*% phiLambda[h:k_ast]
      delta[h] <- rgamma(1, shape=ad, rate=bd)
      # update tau
      tau <- cumprod(delta)
    }
    
    # update Sigma
    res <- (Y - eta %*% t(Lambda))^2
    Sigma_inv <- sapply(1:p, function(j){
      rgamma(1, shape=a_sigma + N / 2, rate=b_sigma + 0.5 * sum(res[,j]))
    })
    Sigma <- diag(1 / Sigma_inv)
    
    ### adaptation strategy ###
    #prob <- exp(alpha0 + alpha1 * iter)
    #u <- runif(1)
    # proportion of elements closed to 0 in each column
    #prop_zero <- apply(Lambda, 2, function(col){
    #  sum(abs(col) < epsilon)/p
    #})
    # index of redundant columns whose elements are all closed to 0
    #idx_red <- which(prop_zero == 1)
    #m <- length(idx_red) # number of redundant factors in each iteration
    
    #if (iter > num_burn && iter %% thin == 0){
    #  k_prev <- dim(cube_eta)[3]
    #  k_curr <- ncol(eta)
    #  if (k_prev < k_curr) {
    #    cube_eta <- abind::abind(cube_eta, 
    #                             array(0, dim=c(N, k_curr-k_prev, num_slice)),
    #                             along=2)
    #    cube_eta[ , , (iter-num_burn)/thin] <- eta
    #  } else {
    #   cube_eta[ , 1:k_curr, (iter-num_burn)/thin] <- eta
    #  }
    #}
    
    #if (u < prob){
    #  if (iter > 20 && m == 0 && all(prop_zero < .995)){
    #    k_ast <- k_ast + 1
    #    Lambda <- cbind(Lambda, rep(0,p))
    #    Z <- cbind(Z, rbinom(N, 1, p0))
    #    eta <- cbind(eta, eta_init(Z[, k_ast], MH, pmom, psi=psi))
    #    phi <- cbind(phi, rgamma(p, shape=df/2, rate=df/2))
    #    delta <- c(delta, rgamma(1, shape=a2, rate=b2))
    #    tau <- cumprod(delta)
    #  } else if (m > 0) {
    #    k_ast <- max(k_ast-m, 1)
    #    Lambda <- Lambda[, -idx_red]
    #    eta <- eta[, -idx_red]
    #    Z <- Z[, -idx_red]
    #    phi <- phi[, -idx_red]
    #    delta <- delta[-idx_red]
    #    tau <- cumprod(delta)
    #  }
    #}
    
    #k_est[iter] <- k_ast
    
    RV_Z[iter] <- RVadj(Z, data[[1]]$Z)
    RV_L[iter] <- RVadj(Lambda, data[[1]]$Lambda)
    RV_L_2[iter] <- RVadj(Lambda%*%t(Lambda), data[[1]]$Lambda%*%t(data[[1]]$Lambda))
    RV_eta[iter] <- RVadj(eta, data[[1]]$eta)
    RV_eta_2[iter] <- RVadj(eta%*%t(eta), data[[1]]$eta%*%t(data[[1]]$eta))
    RV_L_eta[iter] <- RVadj(Lambda%*%t(eta), data[[1]]$Lambda%*%t(data[[1]]$eta))
    
    if (iter %% 100 == 0){
      cat(iter)
      cat("\n")
    }
  }
  
  # Summarize factor score
  sparsity <- apply(cube_eta!=0, 1, function(arr){
    apply(arr, 1, get_mode)
  }) %>% t
  idx_slab <- which(sparsity==1, arr.ind=T)
  idx_spike <- which(sparsity==0, arr.ind=T)
  eta_slab <- sapply(1:nrow(idx_slab), function(m){
    i <- idx_slab[m, 1]
    k <- idx_slab[m, 2]
    mean(cube_eta[i, k, ])
  })
  eta_hat <- matrix(0, N, ncol(sparsity))
  eta_hat[idx_slab] <- eta_slab
  eta_hat[idx_spike] <- 0
  eta_hat <- trim(eta_hat)
  
  #if (is.null(eta0)){
  #  return(list(Lambda=Lambda_hat, eta=eta_hat, k=k_hat, cube_eta=cube_eta))
  #} else {
    # RV <- coeffRV(Lambda_hat, Lambda0)
    # RV <- coeffRV(eta_hat, eta0)
    return(list(RV=RV, eta=eta_hat, cube_eta=cube_eta))
  return(list(RV_Z=RV_Z, RV_L=RV_L, RV_eta=RV_eta,RV_L_2=RV_L_2, RV_eta_2=RV_eta_2, RV_L_eta=RV_L_eta))
  #}
}



plot(RV_Z[100:500], type="l")
plot(RV_L[100:500], type="l")
plot(RV_eta[100:500], type="l")
plot(RV_L_eta[100:500], type="l")
