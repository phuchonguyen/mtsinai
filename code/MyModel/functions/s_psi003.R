s_psi003 <- function(prm, cst, args=NULL) {
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
    osX <- unlist(lapply(1:n, function(i) os[[i]]%*%Xtilda[[i]]))
    oso <- unlist(lapply(1:n, function(i) os[[i]]%*%omega[[i]][,k]))
    Stilda <- rcppeigen_invert_matrix(diag(oso) + Kinv)
    mutilda <- Stilda%*%osX
    psi[,k] <- rcpp_rmvnorm(1, S=Stilda, mu=mutilda)
  }
  
  # Sample eta here instead, so eta doesn't depend on Y
  nu <- array(NA, dim=c(n,K))
  for (i in 1:n) {
    os <- t(omega[[i]])%*%diag(1/prm$sigmax)
    Stilda <- rcppeigen_invert_matrix(diag(1,K) + os%*%omega[[i]])
    mutilda <- Stilda%*%os%*%mu[[i]]
    nu[i,] <- rcpp_rmvnorm(1, S=Stilda, mu=mutilda)
  }
  eta <- psi + nu
  
  prm[["psi"]] <- psi
  prm[["eta"]] <- eta
  return(prm)
}


