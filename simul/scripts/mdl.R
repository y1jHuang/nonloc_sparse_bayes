BFMAN <- function(Y, num_iter = 500,  num_burn = 300, thin = 1, 
                  num_slice = (num_iter - num_burn) / thin, 
                  eta0, Lambda0, Sigma0){
  # Parameter initialization
  N <- nrow(Y)
  p <- ncol(Y)
  #std <- apply(Y, 2, sd)
  # Y <- scale(Y)
  
  k_ast <- ncol(Lambda0) * 2 # Initial guess of k
  psi <- 1 # dispersion parameter for pMOM density
  theta <- rep(0.5, k_ast)
  a_theta <- 1
  b_theta <- 5 # hyperparameter for the prior of theta
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
  
  Z <- matrix(1, N, k_ast)
  eta <- matrix(0, nrow=N, ncol=k_ast)
  Lambda <- matrix(0, nrow=p, ncol=k_ast)
  Sigma_inv <- rgamma(p, shape=a_sigma, rate=b_sigma)
  Sigma <- diag(1/Sigma_inv)
  phi <- rgamma(p*k_ast, shape=df/2, rate=df/2) %>% matrix(nrow=p, ncol=k_ast)
  delta <- c(rgamma(1, shape=a1, rate=b1), rgamma(k_ast-1, shape=a2, rate=b2))
  tau <- cumprod(delta)
  
  # A cube to store eta in each iteration after the burn-in phase
  cube_eta <- array(0, dim=c(N, k_ast, num_slice))
  # A cube to store Lambda in each iteration after the burn-in phase
  cube_Lambda <- array(0, dim=c(p, k_ast, num_slice))
  # A matrix to store averaged eta %*% t(eta)
  Sigma_eta <- matrix(0, N, N)
  # A matrix to store averaged Lambda %*% t(Lambda)
  Sigma_Lambda <- matrix(0, p, p)
  
  # Gibbs sampling
  for (iter in 1:num_iter){
    
    # Update Z
    for (i in 1:N){
      for (k in 1:k_ast){
        a <- 1 / psi + t(Lambda[, k]) %*% diag(Sigma_inv) %*% Lambda[, k]
        b <- t(Y[i,] - Lambda[, -k] %*% eta[i, -k]) %*% diag(Sigma_inv) %*% Lambda[, k]
        t <- (a * psi)^(-3/2) * (1 + b^2 / a) * exp(b^2 / (2 * a))
        # prob <- 1 - (1 - theta) / (1 - theta + theta * t)
        prob <- 1 - (1 - theta[k]) / (1 - theta[k] + theta[k] * t)
        Z[i, k] <- rbinom(1, 1, prob)
      }
    }
    
    # Update eta
    for(i in 1:N){
      for(k in 1:k_ast){
        if (Z[i, k]){
          a <- 0.5 * (Sigma_inv %*% Lambda[,k]^2 + 1/psi)
          b <- 0.5 * Lambda[, k] %*% diag(Sigma_inv) %*% (Lambda[, -k] %*% eta[i,-k] - Y[i,])
          eta[i, k] <- MH(func=pos_eta, init=eta[i, k], M=2, a=a, b=b, sd_value=1)
        } else {
          eta[i, k] <- 0
        }
      }
    }
    
    # Update theta
    theta <- sapply(1:k_ast, function(k){rbeta(1, sum(Z[, k])+a_theta, N-sum(Z[, k])+b_theta)})
    
    # update Lambda
    Lambda <- sapply(1:p, function(j){
      D_inv <- diag(phi[j,] * tau)
      # Covariance of lambda posterior
      cov_lambda <- (D_inv + Sigma_inv[j] * t(eta) %*% eta) %>%
        chol %>% chol2inv
      # Mean of lambda posterior
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
    
    if (iter > num_burn && iter %% thin == 0){
      Sigma_eta <- Sigma_eta + eta %*% t(eta) / num_slice
      Sigma_Lambda <- Sigma_Lambda + Lambda %*% t(Lambda) / num_slice
      cube_Lambda[, , (iter-num_burn)/thin] <- Lambda
      cube_eta[, , (iter-num_burn)/thin] <- eta
    }
    
    if (iter %% 1000 == 0){
      cat(iter)
      cat("\n")
    }
  }

  Lambda_hat <- apply(cube_Lambda, 1, function(arr){
      apply(arr, 1, median)
  }) %>% t
  
  eta_hat <- apply(cube_eta, 1, function(arr){
      apply(arr, 1, median)
  }) %>% t
  
  # Compute RV coefficient
  RV_Lambda = coeffRV(Lambda_hat %*% t(Lambda_hat), Lambda0 %*% t(Lambda0))$rv
  RV_eta = coeffRV(eta_hat %*% t(eta_hat), eta0 %*% t(eta0))$rv
  # RV_Lambda = coeffRV(Lambda_hat, Lambda0)$rv
  # RV_eta = coeffRV(eta_hat, eta0)$rv
  theta_est = (apply(eta_hat!=0, 2, sum) / nrow(eta_hat)) %>%
    sort(decreasing = T)
  return(list(RV_Lambda=RV_Lambda, RV_eta=RV_eta, theta_est=theta_est))
}


MGPS <- function(Y, num_iter = 500,  num_burn = 300, thin = 1,
                 num_slice = (num_iter - num_burn) / thin, 
                 eta0, Lambda0, Sigma0){
  n <- nrow(Y)
  p <- ncol(Y)
  k_tilde <- round(5*log(p)) # initial guess of k
  k_ast <- 10
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
  
  Lambda <- matrix(0, nrow=p, ncol=k_ast)
  sigma_inv <- rgamma(p, shape=a_sigma, rate=b_sigma)
  Sigma <- diag(1/sigma_inv)
  phi <- rgamma(p*k_ast, shape=df/2, rate=df/2) %>% matrix(nrow=p, ncol=k_ast)
  delta <- c(rgamma(1, shape=a1, rate=b1), rgamma(k_ast-1, shape=a2, rate=b2))
  tau <- cumprod(delta)
  
  k_est <- rep(0, num_slice)
  Sigma_eta <- matrix(0, n, n)
  Sigma_Lambda <- matrix(0, p, p)
  
  # Gibbs sampling
  for (iter in 1:num_iter){
    
    # Update eta
      # Covariance of eta posterior
    cov_eta <- (diag(1, k_ast) + t(Lambda) %*% diag(sigma_inv) %*% Lambda) %>%
      chol %>% chol2inv
      # Mean of eta posterior
    mean_eta <- cov_eta %*% t(Lambda) %*% diag(sigma_inv) %*% t(Y)
    eta <- t(mean_eta) + rmvn(n, rep(0, k_ast), cov_eta)
    
    # Update Lambda
    Lambda <- sapply(1:p, function(j){
      D_inv <- diag(phi[j,] * tau)
      # Covariance of lambda posterior
      cov_lambda <- (D_inv + sigma_inv[j] * t(eta) %*% eta) %>%
        chol %>% chol2inv
      # Mean of lambda posterior
      mean_lambda <- (sigma_inv[j] * cov_lambda %*% t(eta) %*% Y[,j]) %>%
        as.vector
      rmvn(1, mean_lambda, cov_lambda)
    }) %>% t
    
    # Update phi
    rate_phi <- df/2 + Lambda^2 %*% diag(tau) # rate of phi posterior
    phi <- sapply(1:p, function(j){
      sapply(1:k_ast, function(h){
        rgamma(1, shape=(df+1)/2, rate=rate_phi[j,h])
      })
    }) %>% t
    
    # Update delta & tau
    phiLambda <- apply(phi*Lambda^2, 2, sum) # the product of phi' and Lambda
    for (h in 1:k_ast){
      if (h == 1){
        ad <- a1 + p * k_ast / 2
      } else {
        ad <- a2 + p * (k_ast - h + 1) / 2
      }
      bd <- 1 + 0.5 / delta[h] * tau[h:k_ast] %*% phiLambda[h:k_ast]
      delta[h] <- rgamma(1, shape=ad, rate=bd)
      tau <- cumprod(delta)
    }
    
    # Update Sigma
    res <- (Y - eta %*% t(Lambda))^2
    sigma_inv <- sapply(1:p, function(j){
      rgamma(1, shape=a_sigma + n / 2, rate=b_sigma + 0.5 * sum(res[,j]))
    })
    Sigma <- diag(1 / sigma_inv)
    
    ### adaptation strategy ###
    prob <- exp(alpha0 + alpha1 * iter)
    u <- runif(1)
    # proportion of elements closed to 0 in each column
    prop_zero <- apply(Lambda, 2, function(col){
      sum(abs(col) < epsilon)/p
    })
    # index of redundant columns whose elements are all closed to 0
    idx_red <- which(prop_zero == 1)
    m <- length(idx_red) # number of redundant factors in each iteration
    
    if (u < prob){
      if (iter > 20 && m == 0 && all(prop_zero < .995)){
        k_ast <- k_ast + 1
        Lambda <- cbind(Lambda, rep(0,p))
        eta <- cbind(eta, rnorm(n, 0, 1))
        phi <- cbind(phi, rgamma(p, shape=df/2, rate=df/2))
        delta <- c(delta, rgamma(1, shape=a2, rate=b2))
        tau <- cumprod(delta)
      } else if (m > 0) {
        k_ast <- max(k_ast-m, 1)
        Lambda <- Lambda[, -idx_red]
        eta <- eta[, -idx_red]
        phi <- phi[, -idx_red]
        delta <- delta[-idx_red]
        tau <- cumprod(delta)
      }
    }
    
    if (iter > num_burn && iter %% thin == 0){
      Sigma_eta <- Sigma_eta + eta %*% t(eta) / num_slice
      Sigma_Lambda <- Sigma_Lambda + Lambda %*% t(Lambda) /num_slice
      mult <- mult + eta %*% t(Lambda) /num_slice
      cov_epsilon <- cov_epsilon + Sigma / num_slice
      k_est[(iter-num_burn)/thin] <- k_ast
    }
    
    if (iter %% 5000 == 0){
      cat(iter)
      cat("\n")
    }
  }
  RV_Lambda <- coeffRV(Sigma_Lambda, Lambda0 %*% t(Lambda0))$rv
  RV_eta <- coeffRV(Sigma_eta, eta0 %*% t(eta0))$rv
  return(list(RV_Lambda=RV_Lambda, RV_eta=RV_eta, k_est=get_mode(k_est)))
}