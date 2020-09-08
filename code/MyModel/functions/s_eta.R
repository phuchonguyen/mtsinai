#TODO: HUGE ISSUE: If Yt's don't have the same dimention, there might be an issue because PcP has different 
#TODO: Recalculate P every time ??!!!
#TODO: There's so much to do to optimize this XD

s_eta <- function(prm, cst, args=NULL) {
  eta <- prm$eta
  
  for (i in unique(cst$idx)) {
    # Propose eta
    #TODO: Code up HMC, for now use MH
    eta_star <- eta
    sigma_eta <- apply(eta, 2, sd)  #TODO: Is there a better way??
    eta_star[i,] <- propose_eta(eta[i,], method="MH", sigma=sigma_eta)
    
    # Accept/Reject proposal  
    eta_loglike <- eta_neg_loglike(eta, prm, cst)
    eta_star_loglike <- eta_neg_loglike(eta_star, prm, cst)
    logr <- eta_loglike - eta_star_loglike
    
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
eta_neg_loglike <- function(eta, prm, cst) {
  etay <- transform_etay(eta, cst$idx, cst$Tx)
  muy <- get_muy(prm, cst, eta=eta)
  
  Y_loglike <- sum(sapply(1:length(muy), 
                          function(i) -0.5*muy[i]^2/prm$sigmay[cst$ty[i]]))
  X_loglike <- sum(sapply(1:nrow(cst$X), 
                          function(i) {
                            mux <- cst$X[i,] - prm$theta%*%prm$xi[i,,]%*%eta[i,]
                            return(-0.5*t(mux)%*%diag(1/prm$sigmax)%*%mux)
                          }))
  eta_loglike <- sum(sapply(1:nrow(eta), 
                            function(i)
                              -0.5*t(eta[i,]-prm$psi[i,])%*%(eta[i,]-prm$psi[i,])))
  return(-(Y_loglike + X_loglike + eta_loglike))
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

