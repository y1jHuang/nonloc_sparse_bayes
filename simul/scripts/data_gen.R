### Data generation
source("../func_lib.R")
# Initialization
num_rep <- 50 # number of replications
N <- 3000 # number of observations
p <- 60 # number of predictors
k <- 6 # true number of factors
theta <- seq(0.8, by = -0.1, length.out = k) # sparsity parameter
num_eff <- seq(2*k, k+1, -1) # number of non-zero elements in each column of Lambda
gen_eta <- function(Z){
  Z <- as.matrix(Z)
  num_slab <- sum(Z)
  eta <- matrix(0, nrow(Z), ncol(Z))
  # rand <- rnorm(num_slab, sd=2)
  rand <- sapply(1:num_slab, function(x){MH(func=pmom, sd_value=0.1, psi=0.5)})
  eta[Z==1] <- rand
  return(eta)
}
data <- lapply(1:num_rep, function(rep){
  Lambda <- LaplacesDemon::rmvn(p, Sigma=diag(k))
  Z <- sapply(theta, function(p){rbinom(N, 1, p)})
  eta <- gen_eta(Z)
  # Sigma <- diag(1/rgamma(p, shape=1, rate=0.25))
  Sigma <- diag(runif(p, 0, 1))
  epsilon <- LaplacesDemon::rmvn(N, 0, Sigma)
  Y <- eta %*% t(Lambda) + epsilon
  return(list(Y=Y, Lambda=Lambda, eta=eta, Sigma=Sigma))
})

# Save data
file_name <- sprintf("../../data/simul/simul_data_p%d_k%d_n%d_rep%d_scenario4.RData", p, k, N, num_rep)
save(data, file=file_name)