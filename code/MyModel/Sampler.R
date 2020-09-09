funcpath <- file.path(getwd(), "code/MyModel/functions")
funcfiles <- list.files(path = funcpath)
sapply(funcfiles, function(x) source(file.path(funcpath, x)))
library("FastGP")
library("Matrix")
library("dplyr")

MySampler <- function(data, niter=5000, nburn=2000, nthin=1,
                      verbose=T) {
  initials <- init_params(data$X, data$Y, 
                        data$ty, data$idy,
                        data$tx, data$idx,
                        data$tau, data$kappa,
                        data$K, data$L)
  params <- initials$params
  const <- initials$const
  
  nout <- floor((niter-nburn)/nthin)
  n <- nrow(const$X)
  p <- ncol(const$X)
  L <- const$L
  K <- const$K
  Tx <- const$Tx
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
    #beta_int = array(NA, dim = c(nout, K*Tx)),
    gamma = array(NA, dim = c(nout, K*Tx)),
    pi_gamma = array(NA, dim = c(nout, 1)),
    sigmay = array(NA, dim = c(nout, Ty)),
    sigmab = array(NA, dim = c(nout, K*Tx))
  )
  
  j <- 1
  at = ceiling(niter/100)
  pb = txtProgressBar(style = 3)
  for (i in 1:niter) {
    # --- verbose --- #
    if (verbose) {
      setTxtProgressBar(pb, i / niter)
    }
    
    # Covariance regression
    params <- s_xi(params, const)
    params <- s_psi(params, const)
    params <- s_theta(params, const)
    params <- s_sigmax(params, const)
    params <- s_phidelta(params, const)
    params <- s_eta(params, const)
    # Y regression
    params <- s_gamma_beta(params, const)
    #TODO remove ?? params[["mu"]] <- get_muy(params, const)
    params <- s_sigmay(params, const)
    
    #TODO: Implement covariates Z
    #params <- update_gamma_beta_z(params, const)
    
    if (i>nburn & i%%nthin==0) {  # TODO: Is this a good way to thin????
      out$eta[j,,] <- params$eta
      out$xi[j,,,] <- params$xi
      out$psi[j,,] <- params$psi
      out$theta[j,,] <- params$theta
      out$sigmax[j,] <- params$sigmax
      out$phi[j,,] <- params$phi
      out$delta[j,] <- params$delta
      #out$beta_int[j,] <- params$beta_int
      out$gamma[j,] <- params$gamma
      out$beta[j,,] <- params$beta
      out$beta0[j,] <- params$beta0
      out$pi_gamma[j,] <- params$pi_gamma
      out$sigmab[j,] <- params$sigmab
      out$sigmay[j,] <- params$sigmay
      
      j <- j + 1
    }
  }
  
  return(out)
}












