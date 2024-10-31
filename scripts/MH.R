################## MH algorithm #####################
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
