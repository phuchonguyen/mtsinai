funcpath <- file.path(getwd(), "code/SelectModel/functions")
funcfiles <- list.files(path = funcpath)
sapply(funcfiles, function(x) source(file.path(funcpath, x)))
library("FastGP")
library("Matrix")
library("dplyr")
library(mvtnorm)
library(fields)

SelectSampler <- function(data, niter=5000, nburn=2000, nthin=1,
                      verbose=T) {
  initials <- init_params(data$X, data$Y, 
                        data$ty, data$idy,
                        data$K)
  prm <- initials$params
  cst <- initials$const
  
  nout <- floor((niter-nburn)/nthin)
  n <- nrow(cst$X)
  p <- ncol(cst$X)
  K <- cst$K
  out <- list(
    beta = array(NA, dim = c(nout, Ty, p)),
    alpha = array(NA, dim = c(nout, Ty)),
    #beta_int = array(NA, dim = c(nout, K*Tx)),
    gamma = array(NA, dim = c(nout, p)),
    pi_gamma = array(NA, dim = c(nout, 1)),
    sigmay = array(NA, dim = c(nout, Ty)),
    sigmab = array(NA, dim = c(nout, p))
  )
  
  j <- 1
  at = ceiling(niter/100)
  pb = txtProgressBar(style = 3)
  for (i in 1:niter) {
    # --- verbose --- #
    if (verbose) {
      setTxtProgressBar(pb, i / niter)
    }
    
    # Y regression
    prm <- s_gamma_beta(prm, cst)  # time elapsed for 10^2 iters: 46.061 *
    prm <- s_sigmay(prm, cst)  # time elapsed for 10^2 iters: 0.780
    
    if (i>nburn & i%%nthin==0) {  # TODO: Is this a good way to thin????
      out$gamma[j,] <- prm$gamma
      out$beta[j,,] <- prm$beta
      out$alpha[j,] <- prm$alpha
      out$pi_gamma[j,] <- prm$pi_gamma
      out$sigmab[j,] <- prm$sigmab
      out$sigmay[j,] <- prm$sigmay
      
      j <- j + 1
    }
  }
  
  return(out)
}












