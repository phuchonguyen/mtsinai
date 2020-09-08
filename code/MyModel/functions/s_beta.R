s_gamma_beta <- function(prm, cst) {
  gamma <- prm$gamma
  beta <- prm$beta
  beta0 <- prm$beta0
  sigmab <- prm$sigmab
  pi_gamma <- prm$pi_gamma
  
  etay <- transform_etay(prm$eta, cst$idx, cst$Tx)
  Syi <- diag(1/prm$sigmay[cst$ty])
  Sy <- diag(prm$sigmay[cst$ty])
  
  p <- length(gamma)
  for (j in 1:p) {
    
    # Sample gamma_j
    muy <- get_muy(prm, cst, j=j)  # Y - eta_-j %*% beta_-j
    Aj <- bdiag(lapply(1:cst$Ty, function(x) etay[,j]))
    betaj <- beta[,j]
    Vg <- rcppeigen_invert_matrix(as.matrix(cst$Vinv/sigmab[j] + t(Aj)%*%Syi%*%Aj))
    Vgd <- det(Vg)
    mj <- as.vector(t(Aj)%*%Syi%*%muy)
    rj <- 1 / exp(-0.5*log(cst$Vdet) + 0.5*log(Vgd) - 0.5*cst$Ty*log(sigmab[j]) + 
             0.5*t(mj)%*%Vg%*%mj) 
    pj <- 1/(1 + (1-pi_gamma)*rj/pi_gamma)
    gamma[j] <- 1*(runif(1) <= pj)  
    
    if (gamma[j] == 1) {
      # Sample beta_j
      mbj <- mj + cst$Vinv%*%rep(1, ncol(cst$Vinv))*beta0[j]/sigmab[j]
      beta[,j] <- as.vector(mvtnorm::rmvnorm(1, mean = Vg%*%mbj, sigma = Vg))
      
      # Sample sigmab_j
      aj <- cst$Ty + cst$asb 
      bj <- t(betaj - beta0[j]) %*% cst$Vinv %*% (betaj - beta0[j]) + cst$bsb
      sigmab[j] <- 1/rgamma(1, 0.5*aj, 0.5*bj)
      
      # Sample beta0_j
      Soj <- rcppeigen_invert_matrix(as.matrix(Aj%*%cst$V%*%t(Aj)*sigmab[j] + Sy))
      IA <- as.vector(Aj%*%rep(1, ncol(Aj)))
      Vo <- 1/as.numeric(t(IA)%*%Soj%*%IA + 1/cst$sb0[j])
      mo <- as.numeric(Vo * t(IA)%*%Soj%*%muy)
      beta0[j] <- rnorm(1, mo, sqrt(Vo))
    }
  }
  
  pi_gamma <- sum(gamma)/length(gamma)
  
  prm[["gamma"]] <- gamma
  prm[["beta"]] <- beta
  prm[["beta0"]] <- beta0
  #prm[["beta_int"]] <- beta_int
  prm[["pi_gamma"]] <- pi_gamma
  prm[["sigmab"]] <- sigmab
  return(prm)
}

