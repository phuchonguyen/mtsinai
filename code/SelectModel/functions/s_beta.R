s_gamma_beta <- function(prm, cst) {
  gamma <- prm$gamma
  beta <- prm$beta
  sigmab <- prm$sigmab
  
  etay <- cst$X
  Syi <- diag(1/prm$sigmay[cst$ty])
  Sy <- diag(prm$sigmay[cst$ty])
  
  p <- length(gamma)
  for (j in 1:p) {
    # Sample gamma_j
    muy <- get_muy_X(prm, cst, j=j)  # Y - eta_-j %*% beta_-j - alpha
    Aj <- diag(1,Ty) %x% etay[,j] 
    Vb <- rcppeigen_invert_matrix(cst$Vinv/sigmab[j] + t(Aj)%*%Syi%*%Aj)
    Vbd <- det(Vb)
    mj <- as.vector(t(Aj)%*%Syi%*%muy)
    rj1 <- -0.5*log(cst$Vdet) + 0.5*log(Vbd) - 0.5*cst$Ty*log(sigmab[j]) +
            0.5*t(mj)%*%Vb%*%mj
    rj1 <- as.numeric(rj1)
    pj0 <- cst$pi_gamma0/(cst$pi_gamma0 + (1-cst$pi_gamma0)*exp(rj1))
    gamma[j] <- 1*(runif(1) > pj0)  
    
    if (gamma[j] == 1) {
      # Sample beta0_j
      # Vg <- Sy + Aj%*%cst$V%*%t(Aj)*sigmab[j]
      # Vgi <- rcppeigen_invert_matrix(Vg)
      # IA <- Aj%*%rep(1, ncol(Aj))
      # Vo <- 1/as.numeric(t(IA)%*%Vgi%*%IA + 1/cst$sb0[j])
      # mo <- as.numeric(Vo * t(IA)%*%Vgi%*%muy)
      # beta0[j] <- rnorm(1, mo, sqrt(Vo))
      
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
  
  prm[["gamma"]] <- gamma
  prm[["beta"]] <- beta
  #prm[["beta_int"]] <- beta_int
  prm[["pi_gamma"]] <- pj0  # Probability that it's 0
  prm[["sigmab"]] <- sigmab
  return(prm)
}

