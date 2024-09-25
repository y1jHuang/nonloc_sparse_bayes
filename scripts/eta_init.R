eta_init <- function(Z, sampler, dens, ...){
  Z <- as.matrix(Z)
  num_slab <- sum(Z)
  eta <- matrix(0, nrow(Z), ncol(Z))
  rand <- sapply(1:num_slab, 
                 function(i){sampler(func=dens, ...)}
                )
  eta[Z==1] <- rand
  return(eta)
}