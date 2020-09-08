s_theta <- function(params, const, args=NULL) {
  X <- const$X
  phi <- params$phi
  delta <- params$delta
  eta <- params$eta
  xi <- params$xi
  theta <- params$theta
  sigmax <- params$sigmax
  
  n <- nrow(X)
  p <- ncol(X)
  tau <- get_factor_tau(delta)
  etatilda <- sapply(1:n, function(i) xi[i,,]%*%eta[i,])
  for (j in 1:p) {
    Stilda <- 1/sigmax[j]*etatilda%*%t(etatilda) + diag(phi[j,]*tau)
    Stilda <- rcppeigen_invert_matrix(Stilda)
    mutilda <- 1/sigmax[j]*Stilda%*%etatilda%*%X[,j]
    theta[j,] <- rcpp_rmvnorm(1, S=Stilda, mu=mutilda)
  }
  
  params[["theta"]] <- theta
  return(params)
}