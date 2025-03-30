library(tidyverse)
library(ggplot2)
library(LaplacesDemon)
source("func_lib.R")

covar <- read.csv("../data/part_derv_lad1.csv") %>%
  filter(MED_ANTIHYPERT == 1 | MED_ANTIDIAB == 1) %>%
  select(PID)
book <- readxl::read_excel("../data/Nutrients_diff.xlsx", sheet = 3)
df <- read.csv("../data/dtia.csv") %>%
  inner_join(covar, by = "PID") %>%
  mutate(SCSFA = DTIA69,
         MCSFA = rowSums(across(DTIA70:DTIA73), na.rm = TRUE),
         LCSFA = rowSums(across(DTIA74:DTIA78), na.rm = TRUE),
         LCMFA = rowSums(across(DTIA80:DTIA83), na.rm = TRUE)) %>%
  select(all_of(book$Variable)) %>%
  filter(if_all(everything(), ~ . >= 0 | is.na(.)))

std <- apply(df, 2, sd) %>% diag
Y <- (df + 0.1) %>%
  log() %>%
  scale(center=T, scale=T)

num_iter=25000
num_burn=10000
thin=1
num_slice=(num_iter - num_burn) / thin

N <- nrow(Y)
p <- ncol(Y)
k_ast <- 12
psi <- 1 # dispersion parameter for pMOM density
theta <- rep(0.3, k_ast)
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

cube_eta <- array(0, dim=c(N, k_ast, num_slice))
cube_Lambda <- array(0, dim=c(p, k_ast, num_slice))

# Gibbs sampling
for (iter in 1:num_iter){
  # Update Z
  for (i in 1:N){
    for (k in 1:k_ast){
      a <- 1 / psi + t(Lambda[, k]) %*% diag(Sigma_inv) %*% Lambda[, k]
      b <- t(Y[i,] - Lambda[, -k] %*% eta[i, -k]) %*% diag(Sigma_inv) %*% Lambda[, k]
      t <- (a * psi)^(-3/2) * (1 + b^2 / a) * exp(b^2 / (2 * a))
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
  # eta <- eta / apply(eta, 2, std)
  # Lambda <- Lambda * apply(eta, 2, std)
  
  # Update theta
  theta <- sapply(1:k_ast, function(k){rbeta(1, sum(Z[, k])+a_theta, N-sum(Z[, k])+b_theta)})
  
  # Update Lambda
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
  Sigma_inv <- sapply(1:p, function(j){
    rgamma(1, shape=a_sigma + N / 2, rate=b_sigma + 0.5 * sum(res[,j]))
  })
  Sigma <- diag(1 / Sigma_inv)
  
  if (iter > num_burn && iter %% thin == 0){
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

# Save data
file_name <- sprintf("../results/nutAnaly_trueK.Rdata")
save(eta_hat, Lambda_hat, file=file_name)