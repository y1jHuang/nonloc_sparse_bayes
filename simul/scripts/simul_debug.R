library(tidyverse)
library(LaplacesDemon)
library(FactoMineR)
source("func_lib.R")

# trim <- function(x){
#   x <- as.matrix(x)
#   idx.tail <- apply(x, 2, function(col){any(col!=0)}) %>%
#     which %>% max
#   return(as.matrix(x[,1:idx.tail]))
# }

file_name <- list.files("../data", pattern="data.*p30_k5_n100.*rep50_norm\\.RData$", full.names=T)
load(file_name)
rep <- 1
Y <- data[[1]]$Y
eta0 <- data[[1]]$eta
Lambda0 <- data[[1]]$Lambda 
Sigma0 <- data[[1]]$Sigma

################### Sparse Bayesian Infinite Factor Model ####################
set.seed(12)
# Lambda0 <- Lambda
num_iter=500
num_burn=250
thin=5
num_slice=(num_iter-num_burn)/thin
# recMethod="PCA"
### Parameter initialization ###
N <- nrow(Y)
p <- ncol(Y)
# std <- apply(Y, 2, sd) # std of Y
# Y <- scale(Y)
k_ast <- ncol(Lambda0)

psi <- 1 # dispersion parameter for pMOM density.
p0 <- 0.6
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

### Gibbs sampling ###
for (iter in 1:num_iter){
  # update Z
  for (i in 1:N){
    for (k in 1:k_ast){
        a <- 1 / psi + t(Lambda[, k]) %*% diag(Sigma_inv) %*% Lambda[, k]
        b <- t(Y[i,] - Lambda[, -k] %*% eta[i, -k]) %*% diag(Sigma_inv) %*% Lambda[, k]
        t <- (a * psi)^(-3/2) * (1 + b^2 / a) * exp(b^2 / (2 * a))
        prob <- 1 - (1 - p0) / (1 - p0 + p0 * t)
        Z[i, k] <- rbinom(1, 1, prob)
    }
  }

  # update eta
  for(i in 1:N){
    for(k in 1:k_ast){
      if (Z[i, k]){
        a <- 0.5 * (Sigma_inv %*% Lambda[,k]^2 + 1/psi)
        b <- 0.5 * Lambda[, k] %*% diag(Sigma_inv) %*% (Lambda[, -k] %*% eta[i,-k] - Y[i,])
        eta[i, k] <- MH(func=pos_eta, step_size=0.1, init=eta[i, k], M=2, a=a, b=b)
      } else {
        eta[i, k] <- 0
      }
    }
  }
  
  # update p0
  p0 <- rbeta(1, sum(Z)+a_p0, N*k_ast-sum(Z)+b_p0)
  
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
  
  if (iter > num_burn && iter %% thin == 0){
    cube_Lambda[, , (iter - num_burn) / thin] <- Lambda
    cube_eta[, , (iter - num_burn) / thin] <- eta
  }
  
  if (iter %% 100 == 0){
    cat(iter)
    cat("\n")
  }
}

RV_Lambda <- coeffRV(cov_Lambda, Lambda0 %*% t(Lambda0))$rv
RV_eta <- coeffRV(cov_eta, eta0 %*% t(eta0))$rv