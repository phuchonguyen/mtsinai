s_sigmax <- function(params, const, args=NULL) {
  X <- const$X
  theta <- params$theta
  xi <- params$xi
  eta <- params$eta
  a0 <- const$asx
  b0 <- const$bsx
  
  p <- ncol(X)
  n <- nrow(X)
  sigmax <- rep(NA, p)
  for (j in 1:p) {
    a <- a0 + n/2
    b <- b0 + sum(sapply(1:n, function(i) (X[i,j] - t(theta[j,])%*%xi[i,,]%*%eta[i,])^2))/2
    sigmax[j] <- 1/rgamma(1, a, b)
  }
  
  params[["sigmax"]] <- sigmax
  return(params)
}