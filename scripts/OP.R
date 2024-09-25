######### OP solution for recovering loading matrix #########
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