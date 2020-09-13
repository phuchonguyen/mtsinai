s_psi <- function(prm, cst, args=NULL) {
  K <- cst$K
  X <- cst$X
  sx <- diag(prm$sigmax)  # sx is a matrix
  psi <- prm$psi
  Kinv <- prm$Kinv
  theta <- prm$theta
  xi <- prm$xi
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
  
  prm[["psi"]] <- psi
  return(prm)
}


