s_theta <- function(prm, cst, args=NULL) {
  
  prm <- s_phidelta(prm, cst)  # time elapsed for 10^2 iters: 0.244
  
  X <- cst$X
  phi <- prm$phi
  delta <- prm$delta
  eta <- prm$eta
  xi <- prm$xi
  theta <- prm$theta
  sigmax <- prm$sigmax
  
  n <- nrow(X)
  p <- ncol(X)
  tau <- get_factor_tau(delta)
  # TODO: quick fix by converting eta[i,] to matrix(eta[i,],...)
  etatilda <- t(sapply(1:n, function(i) xi[i,,]%*% matrix(eta[i,], ncol(eta), 1)))
  ete <- t(etatilda)%*%etatilda
  for (j in 1:p) {
    Stilda <- 1/sigmax[j]*ete + diag(phi[j,]*tau)
    Stilda <- FastGP::rcppeigen_invert_matrix(Stilda)
    mutilda <- 1/sigmax[j]*Stilda%*%t(etatilda)%*%X[,j]
    theta[j,] <- FastGP::rcpp_rmvnorm(1, S=Stilda, mu=mutilda)
  }
  
  prm[["theta"]] <- theta
  return(prm)
}


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
  
  # delta_1 has a different posterior
  tau_1 <- get_factor_tau(delta, exclude_id = 1) 
  b <- sum(sapply(1:L, function(l) tau_1[l]*t(phi[,l])%*%(theta[,l]^2)))
  delta[1] <- rgamma(1, a1+(p*L)/2, 1+0.5*b)
  for (h in 2:L) {
    tau_h <- get_factor_tau(delta, exclude_id = h)
    b <- sum(sapply(1:L, function(l) tau_h[l]*t(phi[,l])%*%(theta[,l]^2)))
    delta[h] <- rgamma(1, a2+p*(L-h+1)/2, 1+0.5*b)
  }
  
  prm[["delta"]] <- delta
  return(prm)
}


s_phi <- function(prm, cst, args=NULL) {
  L <- cst$L
  a0 <- cst$aphi
  b0 <- cst$bphi
  delta <- prm$delta
  theta <- prm$theta
  phi <- prm$phi
  
  p <- nrow(theta)
  tau <- get_factor_tau(delta)
  for (l in 1:L) {
    for (j in 1:p) {
      b <- tau[l]*theta[j,l]^2
      phi[j,l] <- rgamma(1, (1+a0[j,l])/2, (b+b0[j,l])/2)
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