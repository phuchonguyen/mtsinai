s_psi <- function(params, const, args=NULL) {
  K <- const$K
  X <- const$X
  sx <- diag(params$sigmax)  # sx is a matrix
  psi <- params$psi
  Kinv <- params$Kinv
  theta <- params$theta
  xi <- params$xi
  n <- nrow(X)
  
  omega <- lapply(1:n, function(i) theta%*%xi[i,,])
  sinv <- lapply(1:n, function(i) rcppeigen_invert_matrix(omega[[i]]%*%t(omega[[i]]) + sx))
  mu <- lapply(1:n, function(i)  X[i,] - (omega[[i]]%*%psi[i,]))
  
  for (k in 1:K) {
    Xtilda <- lapply(1:n, function(i) mu[[i]] + omega[[i]][,k]*psi[i,k])
    os <- lapply(1:n, function(i) t(omega[[i]][,k])%*%sinv[[i]])
    osX <- lapply(1:n, function(i) os[[i]]%*%Xtilda[[i]])
    oso <- lapply(1:n, function(i) os[[i]]%*%omega[[i]][,k])
    Stilda <- rcppeigen_invert_matrix(as.matrix(bdiag(oso)) + Kinv)
    mutilda <- Stilda%*%unlist(osX)
    psi[,k] <- rcpp_rmvnorm(1, S=Stilda, mu=mutilda)
  }
  
  params[["psi"]] <- psi
  return(params)
}


