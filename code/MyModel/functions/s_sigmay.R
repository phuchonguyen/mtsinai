s_sigmay = function(prm, cst) {
  sigmay <- prm$sigmay
  alpha <- prm$alpha
  
  for (t in 1:cst$Ty) {
    n <- sum(cst$ty==t)
    Yt <- cst$Y[cst$ty==t]
    idyt <- cst$idy[cst$ty==t]
    etayt <- transform_etay(prm$eta, cst$idx, cst$Tx)[idyt,]
    # TODO matrix(...) is hacky to make %*% work when used with scalar
    muy <- Yt - etayt%*%matrix(prm$beta[t,], ncol(prm$beta), 1)
    
    # Sample intercept alpha ~ N(0,10)
    s <- 1/(n/sigmay[t] + 1/10)   
    m <- (sum(muy)/sigmay[t])*s
    alpha[t] <- rnorm(1, m, sqrt(s))
    
    # Sample variance 
    a <- n + cst$asy
    muy <- muy - alpha[t]
    b <- sum(muy^2) + cst$bsy
    sigmay[t] <- 1/rgamma(1, 0.5*a, 0.5*b)
  }
  
  prm[["sigmay"]] = sigmay
  prm[["alpha"]] = alpha
  return(prm)
}
