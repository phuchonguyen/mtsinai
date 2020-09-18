s_eta <- function(prm, cst, args=NULL) {
  eta <- prm$eta
  
  for (i in 1:nrow(eta)) {
    # Propose eta
    #TODO: Code up HMC, for now use MH
    eta_star <- eta
    #eta_star[i,] <- rcpp_rmvnorm(1, mu=prm$psi[i,], S=diag(1, cst$K))  # proposes from the sampling model of X
    eta_star[i,] <- rcpp_rmvnorm(1, mu=eta[i,], S=diag(1, cst$K))
    
    # Accept/Reject proposal  
    eta_l <- eta_loglike(i, eta, prm, cst)
    eta_star_l <- eta_loglike(i, eta_star, prm, cst)
    logr <- eta_star_l - eta_l
    
    if (log(runif(1)) < logr) {
      eta <- eta_star
      prm[["eta_naccept"]] <- prm[["eta_naccept"]] + 1
    }
    
    prm[["eta_npro"]] <- prm[["eta_npro"]] + 1
  }
    
  prm[["eta"]] <- eta
  return(prm)
}


#TODO: Add other covariates
#TODO: Add quadratic terms
eta_loglike <- function(i, eta, prm, cst) {
  id <- cst$idx[i]
  
  etay <- transform_etay(eta, cst$idx, cst$Tx)[id,]
  idc <- cst$idy==id
  muy <- cst$Y[idc] - (diag(1, sum(idc)) %x% t(etay)) %*% as.vector(t(prm$beta)[,cst$ty[idc]]) - prm$alpha[cst$ty[idc]]
  Syi <- diag(1/prm$sigmay[cst$ty[idc]])
  Y_loglike <- -0.5*t(muy)%*%Syi%*%muy
  
  mux <- cst$X[i,] - prm$theta%*%prm$xi[i,,]%*%eta[i,]
  Sxi <- diag(1/prm$sigmax)
  X_loglike <- -0.5*t(mux)%*%Sxi%*%mux
  
  #eta_loglike <- 0 # Since proposes from this
  eta_loglike <- -0.5*t(prm$psi[i,]-eta[i,])%*%(prm$psi[i,]-eta[i,])
  
  return((Y_loglike + X_loglike + eta_loglike))
}

