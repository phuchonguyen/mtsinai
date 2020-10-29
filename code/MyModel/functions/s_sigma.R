s_Sigma <- function(prm, cst) {
  v <- cst$vS
  q <- cst$q
  prm <- s_l(prm, cst, v, q)
  
  Y <- cst$Y
  etay <- transform_etay(prm$eta, cst$idx, cst$Tx) # cst$X  
  n <- cst$N
  p <- sum(prm$gamma)
  B <- prm$B
  l <- prm$l
  
  M <- Y - etay%*%B
  D <- 2*v*diag(1/l)
  #k <- var(as.vector(M))
  #D <- k*diag(q)
  S <- t(M)%*%M + D + t(B)%*%diag(1/prm$nu)%*%B
  nu <- n+p+(v+q-1)
  
  Sigma <- CholWishart::rInvWishart(1, df = nu, Sigma = S)[,,1]
  
  stopifnot(dim(Sigma) == c(cst$q, cst$q))
  
  prm[["Sigma"]] <- Sigma
  prm[["Sigmainv"]] <- chol2inv(chol(Sigma))
  return(prm)
}


s_l <- function(prm, cst, v, q) {
  d <- cst$dS
  Sinv <- prm$Sigmainv
  l <- rep(NA, q)
  for (i in sample(1:q)) {
    l[i] <- 1/rgamma(1, 0.5*(v + q), v*Sinv[i,i] + 1/d[i]^2)
  }
  
  stopifnot(sum(is.na(l))==0)
  
  prm[["l"]] <- l
  return(prm)
}
