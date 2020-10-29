# Adding in Spike and Slab priors for each row of B
# TODO: add option to select groups by factor or exposure time
s_Bsns <- function(prm, cst, prior="laplace") {
  Y    <- cst$Y
  etay <- transform_etay(prm$eta, cst$idx, cst$Tx) # cst$X  
  q    <- cst$q
  B    <- prm$B
  nu   <- prm$nu
  zeta <- prm$zeta
  gamma <- prm$gamma
  pi_gamma <- prm$pi_gamma
  Sinv  <- prm$Sigmainv
  mg    <- 1 # TODO: changes if select groups of variables
  
  for (i in sample(1:nrow(B))) {
    ### Probability of gamma_j = 0 ###
    Xg <- etay[,i]
    Sginv <- as.matrix(1/nu[i] + t(Xg)%*%Xg)
    Sg <- 1/Sginv
    Mg <- Sg %*% t(Xg) %*% (Y - etay[,-i] %*% B[-i,])
    lRg <- log(nu[i])*(-q/2) + log(Sg)*(q/2) +
      0.5*sum(diag( Sinv %*% t(Mg) %*% Sginv %*% Mg))
    p0 <- pi_gamma/(pi_gamma + (1-pi_gamma)*exp(lRg))
    
    ### Update B[i,], nu[i], zeta[i] given gamma[i] ###
    if(runif(1) < p0) {
      
      gamma[i] <- 0
      B[i,] <- 0
      if (prior == "laplace") { 
        # B[i,] ~ MLaplace, nu ~ G((q+1) / 2, zeta / 2)
        nu[i] <- rgamma(1, 0.5*(q+1), 0.5*zeta[i])  # MBSGS prior, MLaplace on Bg
        zeta[i] <- rgamma(1, cst$aB, cst$tauB)    
      } else if (prior == "tpb") {
        # B[i,] ~ TPB, nu ~ G(u, zeta)
        nu[i] <- rgamma(1, cst$uB, zeta[i])
        zeta[i] <- rgamma(1, cst$aB, cst$tauB)
      } else stop(paste("Prior", prior, "not implemented"))
    
    } else {
      
      gamma[i] <- 1
      B[i,] <- MASS::mvrnorm(1, mu = Mg, Sigma = Sg[1,1]*prm$Sigma)
      d <- t(B[i,])%*%Sinv%*%B[i,]
      if (prior == "laplace") {
        nu[i] <- 1/statmod::rinvgauss(1, mean = sqrt(zeta[i]/d), shape = zeta[i])      # MBSGS 
        zeta[i] <- rgamma(1, cst$aB, nu[i]*0.5 + cst$tauB) # zeta ~ G(a, tau) prior
      } else if (prior == "tpb") {
        chi <- max(d, .Machine$double.eps) # to prevent chi parameter from collapsing to 0 from MBSP
        psi <- 2*zeta[i]
        lambda <- cst$uB-q*0.5
        nu[i] <- GIGrvg::rgig(1, lambda = lambda, chi = chi, psi = psi)
        zeta[i] <- rgamma(1, cst$aB, nu[i] + cst$tauB)
      }
    }
    
  }
  
  ### Update P(gamma_j=0): pi_gamma ###
  pi_gamma <- rbeta(1, cst$api + sum(gamma==0), cst$bpi + sum(gamma))
  
  prm[["gamma"]] <- gamma
  prm[["pi_gamma"]] <- pi_gamma
  prm[["B"]]     <- B
  prm[["nu"]]    <- nu
  prm[["zeta"]]  <- zeta
  
  return(prm)
}
