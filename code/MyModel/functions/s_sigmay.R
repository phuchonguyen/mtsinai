s_sigmay = function(prm, cst) {
  
  sigmay <- prm$sigmay
  
  for (t in 1:cst$Ty) {
    n <- sum(cst$ty==t)
    a <- n + cst$asy
    Yt <- cst$Y[cst$ty==t]
    idyt <- cst$idy[cst$ty==t]
    etayt <- transform_etay(prm$eta, cst$idx, cst$Tx)[idyt,]
    muy <- sapply(1:length(Yt), function(i) {
      Yt[i] - t(etayt[i,])%*%(prm$beta[t, ])
    }) 
    b <- t(muy)%*%muy + cst$bsy
    sigmay[t] <- 1/rgamma(1, 0.5*a, 0.5*b)
  }
  
  prm[["sigmay"]] = sigmay
  
  return(prm)
}