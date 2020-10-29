# Sample from GP with Gaussian kernel  
# tau: default unit-scale covariance function
# kappa: bandwidth
# mu: mean function
# v: d x 1 vector of time 
# p: 1 is expo cov function, 2 is square expo cov 
# Returns ========================================
# n independent samples
rgp <- function(n, v, mu, kappa, tau=1, p=2, mu_args=NULL) {
  S <- se_cov(v, sigma=tau, kappa=kappa)
  m <- mu(v)
  X <- MASS::mvrnorm(1, m, S)
  # C <- t(chol(S))
  # # Some mu() returns a matrix, depending on mu_args
  # d <- length(v)
  # Z <- array(rnorm(n*d, 0, 1))
  # X <- m + C%*%Z
  return(X)
}

# Mean function that returns 0
mu_zero <- function(v, args=NULL) {
  v <- as.matrix(v)
  return(array(0, dim = dim(v)))
}

# Mean function returns a  constant
mu_cst <- function(v, args) {
  v <- as.matrix(v)
  return(array(args$cst, dim = dim(v)))
}

# Mean function is linear in v
mu_lin <- function(v, args) {
  v <- as.matrix(v)
  return(args$intercept + args$slope*v)
}


# Sine mean function with correlated variables
mu_sine <- function(v, args=NULL) {
  return(sine(v))
}


disjoint <- function(x){
  if(x <= -0.5) return(x)
  else return(x+1)
}


plot_gp_dif <- function(samps, tx, n, p, Tx) {
  nsamp <- dim(samps$theta)[1]
  n <- length(tx)
  p <- dim(truemu)[2]
  pstmu <- apply(samps$Mu, c(2,3), mean)
  pstSigma <- apply(samps$Psi, c(2,3,4), mean)
  spstmu <- array(NA, dim = dim(truemu))
  spstSigma <- array(NA, dim = dim(trueSigma))
  for (t in 1:Tx) {
    spstmu[t,] <- apply(pstmu[tx==t,], 2, mean)
    spstSigma[t,,] <- apply(pstSigma[tx==t,,], c(2,3) , mean)
  }
  
  if (Tx == 1) {
    ylim <- c(min(-0.5, truemu), max(0.5, truemu))
    par(mfrow=c(1,2))
    plot(1:Tx, truemu[,1],
         ylab="true_mu",
         ylim=ylim, xlim=c(1,Tx), col=1)
    for (j in 2:p) {
      lines(1:Tx, truemu[,j], col=j)
      points(1:Tx, truemu[,j], col=j)
    } 
    ylim <- c(min(-0.5, trueSigma), max(0.5, trueSigma))
    plot(1:Tx, trueSigma[,1,1], 
         ylab="true_Sigma",
         ylim=ylim, xlim=c(1,Tx), col=1)
    for (i in 1:p) {
      for (j in i:p) {
        lines(1:Tx, trueSigma[,i,j], col=(i*j + i))
        points(1:Tx, trueSigma[,i,j], col=(i*j + i))
      }
    }
    
    difmu <- spstmu - truemu
    difSigma <- spstSigma - trueSigma
    p <- dim(trueSigma)[2]
    ylim <- c(min(-0.5, truemu), max(0.5, truemu))
    par(mfrow=c(1,2))
    plot(1:Tx, difmu[,1], 
         ylab="post_mu - true_mu",
         ylim=ylim, xlim=c(1,Tx), col=1)
    for (j in 2:p) {
      lines(1:Tx, difmu[,j], col=j)
      points(1:Tx, difmu[,j], col=j)
    } 
    ylim <- c(min(-0.5, trueSigma), max(trueSigma,.5))
    plot(1:Tx, difSigma[,1,1], 
         ylab="post_Sigma - true_Sigma",
         ylim=ylim, xlim=c(1,Tx), col=1)
    for (i in 1:p) {
      for (j in i:p) {
        lines(1:Tx, difSigma[,i,j], col=(i*j + i))
        points(1:Tx, difSigma[,i,j], col=(i*j + i))
      }
    }
  } else {
    ylim <- c( min(-0.5, truemu), max(truemu, 0.5))
    par(mfrow=c(1,2))
    plot(1:Tx, truemu[,1],
         type="l", 
         ylab="true_mu",
         ylim=ylim, xlim=c(1,Tx), col=1,
         main = "Mean and Covariance for X")
    for (j in 2:p) {
      lines(1:Tx, truemu[,j], col=j)
    } 
    ylim <- c(min(-0.5, trueSigma), max(trueSigma, 0.5))
    plot(1:Tx, trueSigma[,1,1], 
         type="l", 
         ylab="true_Sigma",
         ylim=ylim, xlim=c(1,Tx), col=1)
    for (i in 1:p) {
      for (j in i:p) {
        lines(1:Tx, trueSigma[,i,j], col=(i*j + i))
      }
    }
    
    difmu <- spstmu - truemu
    difSigma <- spstSigma - trueSigma
    p <- dim(trueSigma)[2]
    ylim <- c(min(-0.5, truemu), max(truemu, 0.5))
    par(mfrow=c(1,2))
    plot(1:Tx, difmu[,1], 
         type="l", 
         ylab="post_mu - true_mu",
         ylim=ylim, xlim=c(1,Tx), col=1)
    for (j in 2:p) {
      lines(1:Tx, difmu[,j], col=j)
    } 
    ylim <- c(min(-0.5, trueSigma), max(trueSigma, 0.5))
    plot(1:Tx, difSigma[,1,1], 
         type="l", 
         ylab="post_Sigma - true_Sigma",
         ylim=ylim, xlim=c(1,Tx), col=1)
    for (i in 1:p) {
      for (j in i:p) {
        lines(1:Tx, difSigma[,i,j], col=(i*j + i))
      }
    }
  }
}
