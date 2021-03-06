s_eta <- function(prm, cst, mh_delta=0.05) {
  eta <- prm$eta
  acp <- prm$acp
  
  for (i in 1:nrow(eta)) {
    # Propose eta
    #TODO: Code up HMC, for now use MH
    eta_star <- eta
    eta_star[i,] <- FastGP::rcpp_rmvnorm(1, mu=eta[i,], S=diag(mh_delta, cst$K))
    
    # Accept/Reject proposal  
    eta_l <- eta_loglike(i, eta, prm, cst)
    eta_star_l <- eta_loglike(i, eta_star, prm, cst)
    logr <- eta_star_l - eta_l
    
    if (log(runif(1)) < logr) {
      eta <- eta_star
      acp[i] <- acp[i] + 1
    }
  }
    
  prm[["eta"]] <- eta
  prm[["acp"]] <- acp
  return(prm)
}


#TODO: Add other covariates
#TODO: Add quadratic terms
eta_loglike <- function(i, eta, prm, cst) {
  
  # logY
  id <- cst$idx[i]
  etayi <- transform_etay(eta, cst$idx, cst$Tx)[id,]
  muy <- t(prm$B)%*%etayi
  Y_loglike <- -0.5*(t(muy)%*%prm$Sigmainv%*%muy) 
    + t(muy)%*%prm$Sigmainv%*%cst$Y[id,]
  
  # From old model ---------------------------
  # Y_loglike <- lapply(idy, function(j) {
  #   tj <- cst$ty[j]
  #   muy <- t(etayi)%*%prm$beta[tj,]
  #   res <- -0.5*(muy^2 - 2*muy*cst$Y[j])/prm$sigmay[tj]
  #   return(res)
  # })
  # Y_loglike <- sum(unlist(Y_loglike))
  
  # logX and prior
  etai <- eta[i,]
  Lambda <- prm$theta%*%prm$xi[i,,]   #TODO: maybe save this in prm
  LambdaS <- t(Lambda)%*%diag(1/prm$sigmax)
  X_loglike <- -0.5*(t(etai)%*%(LambdaS%*%Lambda + diag(1, cst$K))%*%etai
                     -2*t(etai)%*%(LambdaS%*%cst$X[i,] + prm$psi[i,]))
  
  return((Y_loglike + X_loglike))
}

