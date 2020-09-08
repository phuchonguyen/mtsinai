s_phidelta <- function(params, const, args=NULL) {
  params <- s_phi(params, const, args=args)
  params <- s_delta(params, const, args=args)
  return(params)
} 


s_delta <- function(params, const, args=NULL) {
  a1 <- const$a1delta
  a2 <- const$a2delta
  L <- const$L
  theta <- params$theta
  phi <- params$phi  # (p x L)
  delta <- params$delta
  p <- nrow(theta)
  
  b <- sapply(1:L, function(l) t(phi[,l])%*%theta[,l]^2)
  # delta_1 has a different posterior
  tau_1 <- get_factor_tau(delta, exclude_id = 1) 
  delta[1] <- rgamma(1, a1+(p*L)/2, 1+(b%*%tau_1)/2)
  for (h in 2:L) {
    tau_h <- get_factor_tau(delta, exclude_id = h)
    delta[h] <- rgamma(1, a2+p*(L-h+1)/2, 1+(b%*%tau_h)/2)
  }
  
  params[["delta"]] <- delta
  return(params)
}


s_phi <- function(params, const, args=NULL) {
  L <- const$L
  b0 <- const$bphi
  delta <- params$delta
  theta <- params$theta
  phi <- params$phi
  
  p <- nrow(theta)
  tau <- get_factor_tau(delta)
  for (j in 1:p) {
    for (l in 1:L) {
      b <- tau[l]*theta[j,l]^2
      phi[j,l] <- rgamma(1, (1+b0)/2, (b+b0)/2)
    }
  }
  
  params[["phi"]] <- phi
  return(params)
}


get_factor_tau <- function(delta, exclude_id=NULL) {
  if (!is.null(exclude_id)) {
    delta[exclude_id] <- 1
  }
  tau <- cumprod(delta)
  return(tau)
}

