s_sigmax <- function(prm, cst, args=NULL) {
  X <- cst$X
  theta <- prm$theta
  xi <- prm$xi
  eta <- prm$eta
  a0 <- cst$asx
  b0 <- cst$bsx
  
  p <- ncol(X)
  n <- nrow(X)
  sigmax <- rep(NA, p)
  for (j in 1:p) {
    a <- a0 + n*0.5
    b <- b0 + 0.5*sum(sapply(1:n, function(i) (X[i,j] - t(theta[j,])%*%xi[i,,]%*%eta[i,])^2))  # TODO: Optimize
    sigmax[j] <- 1/rgamma(1, a, b)
  }
  
  prm[["sigmax"]] <- sigmax
  return(prm)
}