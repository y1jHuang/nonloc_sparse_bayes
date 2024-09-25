pos_eta <- function(x, a=1, b=0, Log=T){
  if(Log){
    - a * (x + b/a)^2 + log(x^2)
  } else {
    exp(- a * (x + b/a)^2) * x^2
  }
}