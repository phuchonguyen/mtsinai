s_xi <- function(prm, cst, args=NULL) {
  L <- cst$L
  K <- cst$K
  X <- cst$X
  theta <- prm$theta
  xi <- prm$xi
  eta <- prm$eta
  sigmax <- prm$sigmax
  Kinv <- prm$Kinv
  n <- nrow(X)
  
  mu <- lapply(1:n, function(i) X[i,] - theta%*%xi[i,,]%*%eta[i,])
  
  for (l in 1:L) {
    for (k in 1:K) {
      sn <- sapply(1:n, function(i) eta[i,k]^2 * sum( theta[,l]^2 / sigmax ))
      Stilda <- FastGP::rcppeigen_invert_matrix(Kinv + diag(sn))
      Xtilda <- lapply(1:n, function(i) mu[[i]] + eta[i,k]*xi[i,l,k]*theta[,l])
      mutilda <- sapply(1:n, function(i) eta[i,k] * sum(theta[,l] * Xtilda[[i]] / sigmax))
      xi[,l,k] <- FastGP::rcpp_rmvnorm(n = 1, S = Stilda, mu = Stilda%*%mutilda)
    }
  }
  
  prm[["xi"]] <- xi
  return(prm)
}