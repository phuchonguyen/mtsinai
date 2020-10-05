s_gamma_beta <- function(prm, cst) {
  gamma <- prm$gamma
  beta <- prm$beta
  sigmab <- prm$sigmab
  pi_gamma <- prm$pi_gamma
  etay <- transform_etay(prm$eta, cst$idx, cst$Tx)
  Syi <- diag(1/prm$sigmay[cst$ty])
  Sy <- diag(prm$sigmay[cst$ty])
  #TODO: Maybe save pj0
  
  p <- length(gamma)
  for (j in 1:p) {
    
    # Sample gamma_j
    muy <- get_muy(prm, cst, j=j)  # Y - eta_-j %*% beta_-j - alpha
    Aj <- diag(1,Ty) %x% etay[,j] 
    Vb <- rcppeigen_invert_matrix(cst$Vinv/sigmab[j] + t(Aj)%*%Syi%*%Aj)
    Vbd <- det(Vb)
    mj <- as.vector(t(Aj)%*%Syi%*%muy)
    rj1 <- -0.5*log(cst$Vdet) + 0.5*log(Vbd) - 0.5*cst$Ty*log(sigmab[j]) +
      0.5*t(mj)%*%Vb%*%mj
    rj1 <- as.numeric(rj1)
    pj0 <- pi_gamma/(pi_gamma + (1-pi_gamma)*exp(rj1))
    gamma[j] <- 1*(runif(1) > pj0)
    
    if (gamma[j] == 1) {
      # Sample beta_j
      mbj <- mj  # + cst$Vinv%*%rep(1, ncol(cst$Vinv))*beta0[j]/sigmab[j]
      beta[,j] <- rcpp_rmvnorm(1, mu = Vb%*%mbj, S = Vb)
      
      # Sample sigmab_j
      aj <- cst$Ty + cst$asb 
      bj <- t(beta[,j])%*%cst$Vinv%*%beta[,j] + cst$bsb
      sigmab[j] <- 1/rgamma(1, 0.5*aj, 0.5*bj)
      
    } else {
      beta[,j] <- 0
      sigmab[j] <- 1/rgamma(1, 0.5*cst$asb, 0.5*cst$bsb)
    }
  }
  
  pi_gamma <- rbeta(1, 1+sum(gamma==0), 1+sum(gamma))
  
  prm[["gamma"]] <- gamma
  prm[["beta"]] <- beta
  #prm[["beta_int"]] <- beta_int
  prm[["pi_gamma"]] <- pi_gamma  # Probability that gammaj is 0
  prm[["sigmab"]] <- sigmab
  return(prm)
}


s_NOgamma_beta <- function(prm, cst) {
  gamma <- prm$gamma
  beta <- prm$beta
  sigmab <- prm$sigmab
  pi_gamma <- prm$pi_gamma
  etay <- transform_etay(prm$eta, cst$idx, cst$Tx)
  Syi <- diag(1/prm$sigmay[cst$ty])
  Sy <- diag(prm$sigmay[cst$ty])
  #TODO: Maybe save pj0
  
  p <- length(gamma)
  for (j in 1:p) {
    
    # Sample gamma_j
    muy <- get_muy(prm, cst, j=j)  # Y - eta_-j %*% beta_-j - alpha
    Aj <- diag(1,Ty) %x% etay[,j] 
    Vb <- rcppeigen_invert_matrix(cst$Vinv/sigmab[j] + t(Aj)%*%Syi%*%Aj)
    Vbd <- det(Vb)
    mj <- as.vector(t(Aj)%*%Syi%*%muy)
    rj1 <- -0.5*log(cst$Vdet) + 0.5*log(Vbd) - 0.5*cst$Ty*log(sigmab[j]) +
      0.5*t(mj)%*%Vb%*%mj
    rj1 <- as.numeric(rj1)
    pj0 <- pi_gamma/(pi_gamma + (1-pi_gamma)*exp(rj1))
    gamma[j] <- 1
    
    if (gamma[j] == 1) {
      # Sample beta_j
      mbj <- mj  # + cst$Vinv%*%rep(1, ncol(cst$Vinv))*beta0[j]/sigmab[j]
      beta[,j] <- rcpp_rmvnorm(1, mu = Vb%*%mbj, S = Vb)
      
      # Sample sigmab_j
      aj <- cst$Ty + cst$asb 
      bj <- t(beta[,j])%*%cst$Vinv%*%beta[,j] + cst$bsb
      sigmab[j] <- 1/rgamma(1, 0.5*aj, 0.5*bj)
      
    } else {
      beta[,j] <- 0
      sigmab[j] <- 1/rgamma(1, 0.5*cst$asb, 0.5*cst$bsb)
    }
  }
  
  pi_gamma <- rbeta(1, 1+sum(gamma==0), 1+sum(gamma))
  
  prm[["gamma"]] <- gamma
  prm[["beta"]] <- beta
  #prm[["beta_int"]] <- beta_int
  prm[["pi_gamma"]] <- pi_gamma  # Probability that gammaj is 0
  prm[["sigmab"]] <- sigmab
  return(prm)
}

