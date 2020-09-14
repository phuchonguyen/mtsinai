# Load
source(file.path(getwd(), "code/CovarModel/CovarModel.R"))


# Constants
SEED <- 123
set.seed(SEED)


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
  C <- t(chol(S))
  # Some mu() returns a matrix, depending on mu_args
  d <- length(v)
  m <- mu(v)
  Z <- array(rnorm(n*d, 0, 1), dim = )
  X <- m + C%*%Z
  return(X)
}


# Mean function that returns 0
mu_zero <- function(v, args=NULL) {
  v <- as.matrix(v)
  return(array(0, dim = dim(v)))
}


# Plotting true functions 
plot_truth <- function(kappa, tau, Tx, L, K, P,
                       a1=10, a2=10, gtheta=3, sxa=2, sxb=0.1) {
  # Generate Theta
  truedelta <- rgamma(1, a1, 1)
  truedelta <- c(truedelta, rgamma(L-1, a2, 1))
  truetau <- get_factor_tau(truedelta)
  truephi <- matrix(rgamma(L*P, gtheta/2, gtheta/2), P, L)
  Theta <- array(NA, dim = c(P, L))
  for (j in 1:P) {
    for (l in 1:L) {
      Theta[j,l] <- rnorm(1, 0, 1/sqrt(truephi[j,l] * truetau[l]))
    }
  }
  # Generate Psi
  truepsi <- cbind(rgp(1, 1:Tx, mu_zero, kappa, tau),
                   rgp(1, 1:Tx, mu_zero, kappa, tau))
  # Generate Xi
  truexi <- array(NA, dim = c(Tx, L, K))
  for (k in 1:K) {
    for (l in 1:L) {
      truexi[,l,k] <- rgp(1, 1:Tx, mu_zero, kappa, tau)
    }
  }
  
  # Use Theta, Xi, Psi to make Mu and Sigma
  truemu <- array(NA, dim = c(Tx, P))
  trueSigma <- array(NA, dim = c(Tx, P, P))
  sigma_x0 <- 1/rgamma(P, sxa, sxb) # TODO: Note this!! In Fox 2015 they use rgamma(1, 0.1) but the variance of that is too big
  for (t in 1:Tx) {
    Lambda <- Theta%*%truexi[t,,]
    truemu[t,] <- Lambda%*%truepsi[t,]
    trueSigma[t,,] <- Lambda%*%t(Lambda)  +  sigma_x0
  }
  ylim <- c(min(truemu), max(truemu))
  par(mfrow=c(1,2))
  plot(1:Tx, truemu[,1], type="l", ylab="true_mu",
       ylim=ylim, xlim=c(1,Tx), col=1)
  for (j in 2:P) {
    lines(1:Tx, truemu[,j], col=j)
  } 
  ylim <- c(min(trueSigma), max(trueSigma))
  plot(1:Tx, trueSigma[,1,1], type="l", ylab="true_Sigma",
       ylim=ylim, xlim=c(1,Tx), col=1)
  for (i in 1:P) {
    for (j in i:P) {
      lines(1:Tx, trueSigma[,i,j], col=(i*j + i))
    }
  }
}

plot_truth(kappa=100, tau=1, Tx=10, L=3, K=2, P=10,
           a1=10, a2=10, gtheta=10)


# Generate data
P <- 10
K <- 2
L <- 3
Tx <- 10
M <- 50     # Number of subjects
KAPPA <- 100 # GP bandwidth
EXP_P <- 2  # GP Gaussian kernel
TAU <- 1    # GP sqrt variance
SXA <- 2
SXB <- 0.1
SIGMA_X0 <- diag(1/rgamma(P, SXA, SXB), P, P)  # Error variance of X
tx <- rep(1:Tx, M)  # TODO: Allow exposures to be missing at different times
idx <- rep(1:M, each=Tx)


# Generate Theta
a1 <- a2 <- gtheta <- 10
truedelta <- rgamma(1, a1, 1)
truedelta <- c(truedelta, rgamma(L-1, a2, 1))
truetau <- get_factor_tau(truedelta)
truephi <- matrix(rgamma(L*P, gtheta/2, gtheta/2), P, L)
Theta <- array(NA, dim = c(P, L))
for (j in 1:P) {
  for (l in 1:L) {
    Theta[j,l] <- rnorm(1, 0, 1/sqrt(truephi[j,l] * truetau[l]))
  }
}


# Generate Psi
truepsi <- cbind(rgp(1, 1:Tx, mu_zero, KAPPA, TAU),
                 rgp(1, 1:Tx, mu_zero, KAPPA, TAU))


# Generate Xi
truexi <- array(NA, dim = c(Tx, L, K))
for (k in 1:K) {
  for (l in 1:L) {
    truexi[,l,k] <- rgp(1, 1:Tx, mu_zero, KAPPA, TAU)
  }
}


# Calculate Mu & Sigma
truemu <- array(NA, dim = c(Tx, P))
trueSigma <- array(NA, dim = c(Tx, P, P))
for (t in 1:Tx) {
  Lambda <- Theta%*%truexi[t,,]
  truemu[t,] <- Lambda%*%truepsi[t,]
  trueSigma[t,,] <- Lambda%*%t(Lambda) + SIGMA_X0
}


# Generate X from the truth
X <- t(sapply(tx, function(t)  mvtnorm::rmvnorm(1, mean=truemu[t,], sigma=trueSigma[t,,])))


# Correlation matrix of X at all times
lattice::levelplot(cor(X))


# Covariance regression
data <- list(
  X=X,
  tx=tx, idx=idx,
  tau=TAU, kappa=KAPPA,
  K=K, L=L
)
niter=10000
nburn=5000
nthin=5
print(paste0("Sampling with niter=", niter, "nburn=", nburn, "nthin=", nthin))
samples <- MySampler(data, niter=niter, nburn=nburn, nthin=nthin)
save(niter, truemu, trueSigma,
     KAPPA, TAU, Tx, idx, tx, X, M,
     SIGMA_X0,
     file=file.path(getwd(),"code/samples/covarXtruth_001.RData"))
saveRDS(samples, file=file.path(getwd(), "code/samples/covarXtruth_001.RDS"))
