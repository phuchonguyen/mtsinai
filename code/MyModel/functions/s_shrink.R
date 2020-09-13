s_phidelta <- function(prm, cst, args=NULL) {
  prm <- s_phi(prm, cst, args=args)
  prm <- s_delta(prm, cst, args=args)
  return(prm)
} 


s_delta <- function(prm, cst, args=NULL) {
  a1 <- cst$a1delta
  a2 <- cst$a2delta
  L <- cst$L
  theta <- prm$theta
  phi <- prm$phi  # (p x L)
  delta <- prm$delta
  p <- nrow(theta)
  
  b <- sapply(1:L, function(l) t(phi[,l])%*%theta[,l]^2)
  # delta_1 has a different posterior
  tau_1 <- get_factor_tau(delta, exclude_id = 1) 
  delta[1] <- rgamma(1, a1+(p*L)/2, 1+(b%*%tau_1)/2)
  for (h in 2:L) {
    tau_h <- get_factor_tau(delta, exclude_id = h)
    delta[h] <- rgamma(1, a2+p*(L-h+1)/2, 1+(b%*%tau_h)/2)
  }
  
  prm[["delta"]] <- delta
  return(prm)
}


s_phi <- function(prm, cst, args=NULL) {
  L <- cst$L
  b0 <- cst$bphi
  delta <- prm$delta
  theta <- prm$theta
  phi <- prm$phi
  
  p <- nrow(theta)
  tau <- get_factor_tau(delta)
  for (j in 1:p) {
    for (l in 1:L) {
      b <- tau[l]*theta[j,l]^2
      phi[j,l] <- rgamma(1, (1+b0)/2, (b+b0)/2)
    }
  }
  
  prm[["phi"]] <- phi
  return(prm)
}


get_factor_tau <- function(delta, exclude_id=NULL) {
  if (!is.null(exclude_id)) {
    delta[exclude_id] <- 1
  }
  tau <- cumprod(delta)
  return(tau)
}

