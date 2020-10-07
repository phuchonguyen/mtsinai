s_theta <- function(prm, cst, args=NULL) {
  X <- cst$X
  phi <- prm$phi
  delta <- prm$delta
  eta <- prm$eta
  xi <- prm$xi
  theta <- prm$theta
  sigmax <- prm$sigmax
  
  n <- nrow(X)
  p <- ncol(X)
  tau <- get_factor_tau(delta)
  # TODO: quick fix by converting eta[i,] to matrix(eta[i,],...)
  etatilda <- t(sapply(1:n, function(i) xi[i,,]%*% matrix(eta[i,], ncol(eta), 1)))
  ete <- t(etatilda)%*%etatilda
  for (j in 1:p) {
    Stilda <- 1/sigmax[j]*ete + diag(phi[j,]*tau)
    Stilda <- rcppeigen_invert_matrix(Stilda)
    mutilda <- 1/sigmax[j]*Stilda%*%t(etatilda)%*%X[,j]
    theta[j,] <- rcpp_rmvnorm(1, S=Stilda, mu=mutilda)
  }
  
  prm[["theta"]] <- theta
  return(prm)
}