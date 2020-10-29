s_B <- function(prm, cst) {
  
  prm <- s_zeta_nu(prm, cst)
  
  Y <- cst$Y
  etay <- transform_etay(prm$eta, cst$idx, cst$Tx) # cst$X  
  ete <- t(etay)%*%etay
  Dinv <- diag(1/prm$nu)
  U <- chol2inv(chol(ete + Dinv))
  V <- prm$Sigma
  M <- U%*%t(etay)%*%Y
  B <- MBSP::matrix.normal(M = M, U = U, V = V)
  
  stopifnot(is.array(B))
  stopifnot(dim(B) == c(cst$Kn, cst$q))
  
  prm[["B"]] <- B
  return(prm)
}


s_zeta_nu <- function(prm, cst) {
  tau <- cst$tauB
  u <- cst$uB
  a <- cst$aB
  q <- cst$q
  nu <- prm$nu
  zeta <- prm$zeta
  Sinv <- prm$Sigmainv
  for (i in sample(1:cst$Kn)) {
    # Option 1: B ~ TPN
    zeta[i] <- rgamma(1, a, nu[i] + tau)
    b <- prm$B[i,]
    chi <- max(t(b)%*%Sinv%*%b, .Machine$double.eps) # to prevent chi parameter from collapsing to 0 from MBSP
    psi <- 2*zeta[i]
    lambda <- u-q*0.5
    nu[i] <- GIGrvg::rgig(1, lambda = lambda, chi = chi, psi = psi)
    # Option 2: B ~ MLaplace
    # d <- t(prm$B[i,])%*%Sinv%*%prm$B[i,]
    # nu[i] <- 1/statmod::rinvgauss(1, mean = sqrt(zeta[i]/d), shape = zeta[i])      # MBSGS 
    # zeta[i] <- rgamma(1, cst$aB, nu[i]*0.5 + cst$tauB)
  }
  
  stopifnot(sum(is.na(zeta)) == 0)
  stopifnot(sum(is.na(nu)) == 0)
  
  prm[["nu"]] <- nu
  prm[["zeta"]] <- zeta
  return(prm)
}
