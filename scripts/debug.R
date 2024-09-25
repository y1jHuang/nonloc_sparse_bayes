library(tidyverse)
get_mode <- function(x){
  # function  used for obtaining mode from a series of numbers
  uniq_vals <- unique(x)
  freq <- tabulate(match(x, uniq_vals))
  uniq_vals[which.max(freq)]
}
N <- 100
k <- 5
num_iter <- 1000
Z <- array(0, dim=c(N,k, num_iter))
Z[,1:(k-2),999] <- matrix(NA, N, k-2)
# rand <- array(rbinom(N*k*600,1, 0.6), dim=c(N,k,600))
# Z[,,1:600] <- rand
# Z <- array(rbinom(N*k*num_iter, 1, 0.6), dim=c(N, k, num_iter))
sparsity <- apply(Z, 1, function(arr){
  apply(arr, 1, get_mode)
}) %>% t
a <- matrix(rnorm(N*k), N, k)
a[sample(1:100, 50), ] <- 0
b <- matrix(0, N, 2)
c <- cbind(b, a, b, a, a, b, a, b)
d <- trim(c)
a <- apply(cube_eta!=0, 1, function(arr){
  apply(arr, 1, sum)
}) %>% t
Z <- matrix(1:100, 20, 5)
sapply(2:5, function(i){
  Z[,i] <- Z[,i-1] + Z[,i]
  return()
})

x=seq(0, 0.5, 0.01)
y=dgamma(x, shape=5, scale=10)
rinvgamma(10000, shape=5, scale=10) %>% hist(breaks=100, freq=F)
sample <- (1/rgamma(10000, shape=5, scale=10)) %>% hist(breaks=100, freq=F)
plot(x, y, type='l')

