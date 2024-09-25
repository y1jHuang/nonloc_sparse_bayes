pmom <- function(x, psi=0.01, Log=T){
  if(Log){
    - 0.5 * log(2 * pi) - 1.5 * log(psi) - 0.5 * x^2 / psi + log(x^2)
  } else {
    (2 * pi)^(-0.5) * psi^(-1.5) * exp(-0.5 * x^2 / psi) * x^2
  }
}