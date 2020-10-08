funcpath <- file.path(getwd(), "MyModel/functions")
funcfiles <- list.files(path = funcpath)
sapply(funcfiles, function(x) source(file.path(funcpath, x)))
library(FastGP)
library(Matrix)
library(dplyr)
library(mvtnorm)
library(fields)


# Update MH for eta with adjusted step size
MySampler <- function(data, niter=5000, nburn=2000, nthin=1,
                      verbose=T) {
  initials <- init_params(data$X, data$Y, 
                          data$ty, data$idy,
                          data$tx, data$idx,
                          data$tau, data$kappa,
                          data$K, data$L)
  prm <- initials$params
  cst <- initials$const
  
  nout <- floor((niter-nburn)/nthin)
  n <- nrow(cst$X)
  p <- ncol(cst$X)
  L <- cst$L
  K <- cst$K
  Tx <- cst$Tx
  out <- list(
    eta = array(NA, dim = c(nout, n, K)),
    xi = array(NA, dim = c(nout, n, L, K)),
    psi = array(NA, dim = c(nout, n, K)),
    theta = array(NA, dim = c(nout, p, L)),
    sigmax = array(NA, dim = c(nout, p)),
    phi = array(NA, dim = c(nout, p, L)),
    delta = array(NA, dim = c(nout, L)),
    beta = array(NA, dim = c(nout, Ty, K*Tx)),
    beta0 = array(NA, dim = c(nout, K*Tx)),
    alpha = array(NA, dim = c(nout, Ty)),
    #beta_int = array(NA, dim = c(nout, K*Tx)),
    gamma = array(NA, dim = c(nout, K*Tx)),
    pi_gamma = array(NA, dim = c(nout, 1)),
    sigmay = array(NA, dim = c(nout, Ty)),
    sigmab = array(NA, dim = c(nout, K*Tx))
  )
  
  j <- 1
  at = ceiling(niter/100)
  pb = txtProgressBar(style = 3)
  # TODO: set this mh_delta in init_param
  mh_delta = 0.05
  for (i in 1:niter) {
    # --- verbose --- #
    if (verbose) {
      setTxtProgressBar(pb, i / niter)
    }
    
    # Covariance regression
    prm <- s_xi(prm, cst)  # time elapsed for 10^2 iters: 50.812 *
    prm <- s_psi(prm, cst)  # time elapsed for 10^2 iters: 26.589 *
    prm <- s_theta(prm, cst)  # time elapsed for 10^2 iters: 0.449 
    prm <- s_sigmax(prm, cst)  # time elapsed for 10^2 iters: 5.878
    prm <- s_phidelta(prm, cst)  # time elapsed for 10^2 iters: 0.244
    # TODO: parallelize this eta step
    prm <- s_eta(prm, cst, mh_delta)  # time elapsed for 10^2 iters: ?? 
    # Y regression
    prm <- s_gamma_beta(prm, cst)  # time elapsed for 10^2 iters: 46.061 *
    prm <- s_sigmay(prm, cst)  # time elapsed for 10^2 iters: 0.780
    
    # s <- proc.time()
    # for (i in 1:10^2) {
    # }
    # print(proc.time()-s)
    
    #TODO: Implement covariates Z
    #prm <- update_gamma_beta_z(prm, cst)
    
    # Update step_size for s_eta
    if (i%%100==0 & i<=nburn) {
      acp_mean = mean(prm$acp)/100
      if(acp_mean > 0.3) {
        mh_delta = mh_delta*2
        print(paste("Updated mh_delta: ", mh_delta))
      } else if(acp_mean < 0.2) {
        mh_delta = mh_delta*2/3
        print(paste("Updated mh_delta: ", mh_delta))
        }
      prm[["acp"]] = rep(0, nrow(prm$eta))
      print(paste("Mean acp: ", acp_mean))
    }
    
    
    if (i>nburn & i%%nthin==0) {  # TODO: Is this a good way to thin????
      out$eta[j,,] <- prm$eta
      out$xi[j,,,] <- prm$xi
      out$psi[j,,] <- prm$psi
      out$theta[j,,] <- prm$theta
      out$sigmax[j,] <- prm$sigmax
      out$phi[j,,] <- prm$phi
      out$delta[j,] <- prm$delta
      #out$beta_int[j,] <- prm$beta_int
      out$gamma[j,] <- prm$gamma
      out$beta[j,,] <- prm$beta
      out$beta0[j,] <- prm$beta0
      out$alpha[j,] <- prm$alpha
      out$pi_gamma[j,] <- prm$pi_gamma
      out$sigmab[j,] <- prm$sigmab
      out$sigmay[j,] <- prm$sigmay
      
      j <- j + 1
    }
  }
  
  return(out)
}


# Regrssion on Y only, given eta
MySamplerNoGamma <- function(data, eta, niter=5000, nburn=2000, nthin=1,
                         verbose=T) {
  initials <- init_params(data$X, data$Y, 
                          data$ty, data$idy,
                          data$tx, data$idx,
                          data$tau, data$kappa,
                          data$K, data$L)
  prm <- initials$params
  cst <- initials$const
  
  nout <- floor((niter-nburn)/nthin)
  n <- nrow(cst$X)
  p <- ncol(cst$X)
  L <- cst$L
  K <- cst$K
  Tx <- cst$Tx
  out <- list(
    eta = array(NA, dim = c(nout, n, K)),
    xi = array(NA, dim = c(nout, n, L, K)),
    psi = array(NA, dim = c(nout, n, K)),
    theta = array(NA, dim = c(nout, p, L)),
    sigmax = array(NA, dim = c(nout, p)),
    phi = array(NA, dim = c(nout, p, L)),
    delta = array(NA, dim = c(nout, L)),
    beta = array(NA, dim = c(nout, Ty, K*Tx)),
    beta0 = array(NA, dim = c(nout, K*Tx)),
    alpha = array(NA, dim = c(nout, Ty)),
    #beta_int = array(NA, dim = c(nout, K*Tx)),
    gamma = array(NA, dim = c(nout, K*Tx)),
    pi_gamma = array(NA, dim = c(nout, 1)),
    sigmay = array(NA, dim = c(nout, Ty)),
    sigmab = array(NA, dim = c(nout, K*Tx))
  )
  
  j <- 1
  at = ceiling(niter/100)
  pb = txtProgressBar(style = 3)
  # TODO: set this mh_delta in init_param
  mh_delta = 0.05
  for (i in 1:niter) {
    # --- verbose --- #
    if (verbose) {
      setTxtProgressBar(pb, i / niter)
    }
    
    # Covariance regression
    prm <- s_xi(prm, cst)  # time elapsed for 10^2 iters: 50.812 *
    # NOTE: samples eta in psi without dependence on Y
    prm <- s_psi(prm, cst)  # time elapsed for 10^2 iters: 26.589 *
    prm <- s_theta(prm, cst)  # time elapsed for 10^2 iters: 0.449 
    prm <- s_sigmax(prm, cst)  # time elapsed for 10^2 iters: 5.878
    prm <- s_phidelta(prm, cst)  # time elapsed for 10^2 iters: 0.244
    # TODO: parallelize this eta step
    # Y regression
    prm[["eta"]] <- eta
    prm <- s_NOgamma_beta(prm, cst)  # time elapsed for 10^2 iters: 46.061 *
    prm <- s_sigmay(prm, cst)  # time elapsed for 10^2 iters: 0.780
    
    # s <- proc.time()
    # for (i in 1:10^2) {
    # }
    # print(proc.time()-s)
    
    #TODO: Implement covariates 
    #prm <- update_gamma_beta_z(prm, cst)
    
    if (i>nburn & i%%nthin==0) {  # TODO: Is this a good way to thin????
      out$eta[j,,] <- prm$eta
      out$xi[j,,,] <- prm$xi
      out$psi[j,,] <- prm$psi
      out$theta[j,,] <- prm$theta
      out$sigmax[j,] <- prm$sigmax
      out$phi[j,,] <- prm$phi
      out$delta[j,] <- prm$delta
      #out$beta_int[j,] <- prm$beta_int
      out$gamma[j,] <- prm$gamma
      out$beta[j,,] <- prm$beta
      out$beta0[j,] <- prm$beta0
      out$alpha[j,] <- prm$alpha
      out$pi_gamma[j,] <- prm$pi_gamma
      out$sigmab[j,] <- prm$sigmab
      out$sigmay[j,] <- prm$sigmay
      
      j <- j + 1
    }
  }
  
  return(out)
}













