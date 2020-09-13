#TODO: HUGE ISSUE: If Yt's don't have the same dimention, there might be an issue because PcP has different 
#TODO: Recalculate P every time ??!!!
#TODO: There's so much to do to optimize this XD

s_eta <- function(prm, cst, args=NULL) {
  eta <- prm$eta
  
  for (i in 1:nrow(eta)) {
    # Propose eta
    #TODO: Code up HMC, for now use MH
    eta_star <- eta
    eta_star[i,] <- mvtnorm::rmvnorm(1, mean=prm$psi[i,])  # proposes from the sampling model
    
    # Accept/Reject proposal  
    eta_l <- eta_loglike(i, eta, prm, cst)
    eta_star_l <- eta_loglike(i, eta_star, prm, cst)
    logr <- eta_star_l - eta_l
    
    if (log(runif(1)) < logr) {
      eta <- eta_star
      #TODO: save acceptance prob
    }
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
  
  eta_loglike <- -0.5*t(prm$psi[i,]-eta[i,])%*%(prm$psi[i,]-eta[i,])
  
  return((Y_loglike + X_loglike + eta_loglike))
}


# https://colindcarroll.com/2019/04/11/hamiltonian-monte-carlo-from-scratch/
propose_eta <- function(eta, method="HMC", V=NULL, p=NULL, sigma=NULL) {
  if (method == "HMC") {
    if (is.null(V) | is.null(p)) stop("When method is HMC, V must not be null")
    stop("Not implemented")
    0
  }
  
  if (method == "MH") {
    eta_star <- rnorm(length(eta), eta, sigma)
    return(eta_star)
  }
}

