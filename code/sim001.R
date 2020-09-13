#' ---
#' title: "Simulation Study"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---
print(getwd())
source("/work/phn5/mtsinai/mtsinai/code/MyModel/Sampler.R")
# TODO: By the model notation: N is number of cases, S is N*Tx
# TODO: Use function to_etai_reg in helpers.R 

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

# Plotting true functions 
plot_truth <- function(kappa, tau, Tx, L, K, P) {
  Theta <- array(rnorm(P*L, 0, 1), dim = c(P, L))
  truepsi <- cbind(rgp(1, 1:Tx, mu_zero, kappa, tau),
                   rgp(1, 1:Tx, mu_zero, kappa, tau))
  truexi <- array(NA, dim = c(Tx, L, K))
  for (k in 1:K) {
    for (l in 1:L) {
      truexi[,l,k] <- rgp(1, 1:Tx, mu_zero, kappa, tau)
    }
  }
  truemu <- array(NA, dim = c(Tx, P))
  trueSigma <- array(NA, dim = c(Tx, P, P))
  for (t in 1:Tx) {
    Lambda <- Theta%*%truexi[t,,]
    truemu[t,] <- Lambda%*%truepsi[t,]
    trueSigma[t,,] <- Lambda%*%t(Lambda) + diag(1, P, P)
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

plot_truth(kappa=100, tau=1, Tx=10, L=3, K=2, P=10)

#' # Simulate data from Gaussian Process in mean and variance for X
#' $$X_{it} = \Lambda(t)\eta_{it}$$
#' $$\Lambda(t) = 1 \forall t $$
#' $$\eta_{it} = t \texttt{ or } a_0, a_0 \in R$$
#' n=10000, T=5, T*=10, p=10, k=4, L=5, 
#' $\kappa_\psi=\kappa=2 (bandwidth param wants to be small since the data wiggles fast from t1 to t2)
#' Assume sites are exchangeable for now...
#' Evaluate estimated $\mu(t)= \Lambda(t)\psi(t), \Sigma(t) = \Lambda(t)\Lambda(t)^T + \Sigma_0^X$ 
#' using Frobenius norm.
#+ cache=TRUE
P <- 10
K <- 2
L <- 3
Tx <- 10
M <- 50     # Number of subjects
KAPPA <- 100 # GP bandwidth
EXP_P <- 2  # GP Gaussian kernel
TAU <- 1    # GP sqrt variance
SIGMA_X0 <- diag(1, P, P)  # Error variance of X
Theta <- array(rnorm(P*L, 0, 1), dim = c(P, L))
tx <- rep(1:Tx, M)  # TODO: Allow exposures to be missing at different times
idx <- rep(1:M, each=Tx)
truepsi <- cbind(rgp(1, 1:Tx, mu_zero, KAPPA, TAU),
                 rgp(1, 1:Tx, mu_zero, KAPPA, TAU))
truexi <- array(NA, dim = c(Tx, L, K))
for (k in 1:K) {
  for (l in 1:L) {
    truexi[,l,k] <- rgp(1, 1:Tx, mu_zero, KAPPA, TAU)
  }
}
truemu <- array(NA, dim = c(Tx, P))
trueSigma <- array(NA, dim = c(Tx, P, P))
for (t in 1:Tx) {
  Lambda <- Theta%*%truexi[t,,]
  truemu[t,] <- Lambda%*%truepsi[t,]
  trueSigma[t,,] <- Lambda%*%t(Lambda) + SIGMA_X0
}

#' ### Generate X from the truth

psi <- t(sapply(tx, function(t) truepsi[t,]))
xi <- array(unlist(lapply(tx, function(t) truexi[t,,])), dim=c(L, K, Tx*M))
xi <- aperm(xi, c(3,1,2))
eta <- t(apply(psi, 1, function(x) mvtnorm::rmvnorm(1, mean=x)))
epsilon_x <- mvtnorm::rmvnorm(length(tx), mean = rep(0, P), sigma = SIGMA_X0)
X <- t(sapply(1:length(tx), function(i) Theta%*%xi[i,,]%*%eta[i,] + epsilon_x[i,]))
print("Simulated X")

#' Correlation matrix of X at all times
lattice::levelplot(cor(X))


#' ### Some scenarios for Y with only the first half of etas are included:
#' 1. Linear effects varies by time
#' $$y = \beta \eta$$
#' 2. TODO: Linear effects varies by time and interactions
#' $$y = \beta \eta + \eta^T \Delta \eta$$

# Each row i is set of coefs for Ty=i for all factors at all time Tx
Ty <- 5
ty <- rep(1:Ty, each=M)
idy <- rep(1:M, Ty)
SIGMA_Y <- rgamma(Ty, 1, 1)
GAMMA <- purrr::rbernoulli(K*Tx, p=0.3)*1
BETA0 <- rnorm(K*Tx, 10, 1) * GAMMA
SIGMA_B <- rgamma(K*Tx, 1, 1)
BETA <- array(NA, dim = c(Ty, K*Tx))
BETA[1,] <- rnorm(K*Tx, mean = BETA0, sd = SIGMA_B) * GAMMA
for (t in 2:Ty) {
  BETA[t,] <- rnorm(K*Tx, mean = BETA[t-1,], sd = SIGMA_B) * GAMMA
}
# TODO add Interactions
etay <- transform_etay(eta, idx, Tx)
Y <- rep(NA, Ty*M)
for (t in 1:Ty) {
  Y[((t-1)*M + 1) : (t*M)] <- etay%*%BETA[t,] + rnorm(M, 0, SIGMA_Y[t])
}
par(mfrow=c(3,3))
plot(etay[,1], Y[1:M], col=GAMMA[1]+1, cex=1)
plot(etay[,3], Y[1:M], col=GAMMA[3]+1, cex=1)
plot(etay[,5], Y[1:M], col=GAMMA[5]+1, cex=1)
plot(etay[,2], Y[1:M], col=GAMMA[2]+1, cex=1)
plot(etay[,9], Y[1:M], col=GAMMA[9]+1, cex=1)
plot(etay[,20], Y[1:M], col=GAMMA[20]+1, cex=1)
plot(etay[,2], Y[(1+M):(2*M)], col=GAMMA[2]+1, cex=1)
plot(etay[,9], Y[(1+M):(2*M)], col=GAMMA[9]+1, cex=1)
plot(etay[,20], Y[(1+M):(2*M)], col=GAMMA[20]+1, cex=1)

par(mfrow=c(1,2))
plot(1:Ty, BETA[,1])
abline(h=BETA0[1], col="red")
plot(1:Ty, BETA[,2])
abline(h=BETA0[2], col="red")

print("Simulated Y")

#' ### Running my sampler
data <- list(
  X=X, Y=Y, 
  ty=ty, idy=idy,
  tx=tx, idx=idx,
  tau=TAU, kappa=KAPPA,
  K=K, L=L
)

niter=10000
nburn=5000
nthin=5
print(paste0("Sampling with niter=", niter, "nburn=", nburn, "nthin=", nthin))
samples <- MySampler(data, niter=niter, nburn=nburn, nthin=nthin)
save(niter, truemu, trueSigma, truepsi, truexi, Theta, KAPPA, TAU, Tx, Ty, idx, tx, idy, ty, X, M,
     Y, BETA0, BETA, GAMMA, SIGMA_B, SIGMA_Y, SIGMA_X0, file="/work/phn5/mtsinai/mtsinai/code/samples/lintruth_006.RData")
saveRDS(samples, file="/work/phn5/mtsinai/mtsinai/code/samples/linsamples_006.RDS")
print("Done")
