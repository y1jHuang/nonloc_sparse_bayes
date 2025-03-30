# Density of pMOM distribution
pmom <- function(x, psi=0.01, Log=T){
  if(Log){
    - 0.5 * log(2 * pi) - 1.5 * log(psi) - 0.5 * x^2 / psi + log(x^2)
  } else {
    (2 * pi)^(-0.5) * psi^(-1.5) * exp(-0.5 * x^2 / psi) * x^2
  }
}

# Posterior of eta
pos_eta <- function(x, a=1, b=0, log=T){
  if(log){
    - a * (x + b/a)^2 + log(x^2)
  } else {
    exp(- a * (x + b/a)^2) * x^2
  }
}

# Get mode from a vector
get_mode <- function(x){
  uniq_vals <- unique(x)
  freq <- tabulate(match(x, uniq_vals))
  uniq_vals[which.max(freq)]
}

# Metropolis-Hastings algorithm
MH <- function(init=0, M=200, num_burn=100, step_size=0.1, func, ...){
  param <- c(init, rep(0, M-1))
  for(i in 2:M){
    cand <- rnorm(n=1, mean=param[i-1], sd=step_size)
    r <- exp(func(cand, ...) - func(param[i-1], ...))
    u <- runif(1, 0, 1)
    z <- u <= r
    param[i] <- z * cand + (1-z) * param[i-1]
  }
  return(param[M])
}

# Approach to solve orthogonal procrustes problem
OP <- function(L_list, tol=10^-5, itermax=50){
  iter <- 1
  L_ast <- L_list[[length(L_list)]]
  repeat {
    if (iter > itermax) break
    L_tild <- sapply(L_list, 
                     function(L, L_ast){
                       s <- svd(t(L) %*% L_ast)
                       return(L %*% s$u %*% t(s$v))},
                     L_ast) %>% t %>%
      apply(2, mean) %>% 
      matrix(nrow(L_ast), ncol(L_ast))
    dist <- norm(L_tild - L_ast, type="F")
    iter <- iter + 1
    L_ast <- L_tild
    if (dist <= tol) break
  }
  return(L_tild)
}