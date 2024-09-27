### Data generation
source("eta_init.R")
source("MH.R")
source("pmom.R")
# Initialization
num_rep <- 50 # number of replications
N <- 100 # number of observations
p <- 30 # number of predictors
k <- 5 # true number of factors
p0 <- 0.6 # probability that factor score is non-zero
# num_eff <- k + sample(k) # number of non-zero elements in each column of Lambda
num_eff <- seq(2*k, k+1, -1)
gen_eta <- function(Z){
  Z <- as.matrix(Z)
  num_slab <- sum(Z)
  eta <- matrix(0, nrow(Z), ncol(Z))
  # rand <- rnorm(num_slab)
  rand <- sapply(1:num_slab, function(x){MH(func=pmom, sd_value=0.1, psi=0.5)})
  eta[Z==1] <- rand
  return(eta)
}
rbiunif <- function(n=1, min=0.7, max=1){
  sign <- rbinom(n, 1, 0.5)
  (sign - 0.5) * 2 * runif(n, min, max)
}
data <- lapply(1:num_rep, function(rep){
  # Simulation for Lambda, Sigma and epsilon
  Lambda <- sapply(1:k, function(i){
    lambda <- rep(0, p)
    idx <- sample(1:p, num_eff[i]) # index of non-zero elements in each column
    lambda[idx] <- rbiunif(num_eff[i], 0.7, 1)
    # lambda[idx] <- rnorm(num_eff[i], 0, 9) # assign values to non-zero elements
    return(lambda)
  })
  Z <- matrix(rbinom(N*k, 1, p0), N, k)
  eta <- gen_eta(Z)
  # Sigma <- diag(1/rgamma(p, shape=1, rate=0.25))
  Sigma <- diag(runif(p, 0, 0.5))
  epsilon <- LaplacesDemon::rmvn(N, 0, Sigma)
  Y <- eta %*% t(Lambda) + epsilon
  return(list(Y=Y, Lambda=Lambda, eta=eta))
})

# Save data
file_name <- sprintf("../data/simul_data_p%d_k%d_n%d_rep%d_norm.RData", p, k, N, num_rep)
save(data, file=file_name)
