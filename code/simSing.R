.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))
source(file.path(getwd(), "MyModel/Sampler.R"))

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

#+ cache=TRUE
P <- 10
K <- 1
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
truepsi <- array(NA, dim = c(Tx, K))
for (i in 1:K) {
  truepsi[,i] <- rgp(1, 1:Tx, mu_zero, KAPPA, TAU)
}

# Generate Xi
truexi <- array(NA, dim = c(Tx, L, K))
for (k in 1:K) {
  for (l in 1:L) {
    truexi[,l,k] <- rgp(1, 1:Tx, mu_zero, KAPPA, TAU)
  }
}


# Generate eta
psi <- matrix(sapply(tx, function(t) truepsi[t,]), M*Tx, K)
eta <- matrix(apply(psi, 1, function(x) mvtnorm::rmvnorm(1, mean=x)), M*Tx, K)


# Calculate Mu & Sigma
truemu <- array(NA, dim = c(Tx, P))
trueSigma <- array(NA, dim = c(Tx, P, P))
for (t in 1:Tx) {
  Lambda <- Theta%*%truexi[t,,]
  truemu[t,] <- Lambda%*%truepsi[t,]
  trueSigma[t,,] <- Lambda%*%t(Lambda) + SIGMA_X0
}

#' ### Generate X from the truth
X <- t(sapply(tx, function(t)  mvtnorm::rmvnorm(1, mean=truemu[t,], sigma=trueSigma[t,,])))
print("Simulated X")

# Each row i is set of coefs for Ty=i for all factors at all time Tx
Ty <- 5
ty <- rep(1:Ty, each=M)
idy <- rep(1:M, Ty)
SIGMA_Y <- sort(1/rgamma(Ty, 10, 1)) 
GAMMA <- (runif(K*Tx) <= 0.3)*1
SIGMA_B <- 1/rgamma(K*Tx, 5, 1)
BETA <- array(NA, dim = c(Ty, K*Tx))
BETA[1,] <- rnorm(K*Tx, mean = rep(1, K*Tx), sd = sqrt(SIGMA_B)) * GAMMA
for (t in 2:Ty) {
  BETA[t,] <- rnorm(K*Tx, mean = BETA[t-1,], sd = sqrt(SIGMA_B)) * GAMMA
}
ALPHA <- rnorm(Ty)
# TODO add Interactions
etay <- transform_etay(eta, idx, Tx)
Y <- rep(NA, Ty*M)
for (t in 1:Ty) {
  Y[ty==t] <- ALPHA[t] + etay%*%BETA[t,] + rnorm(sum(ty==t), 0, sqrt(SIGMA_Y[t]))
}

print("Simulated Y")


### Running my sampler for Y + X
data <- list(
  X=X, Y=Y,
  ty=ty, idy=idy,
  tx=tx, idx=idx,
  tau=TAU, kappa=KAPPA,
  K=K, L=L
)

niter=20000
nburn=10000
nthin=5
print(paste0("Sampling with niter = ", niter, " nburn = ", nburn, " nthin = ", nthin))
samples <- MySampler003(data, niter=niter, nburn=nburn, nthin=nthin)
filename <- "lintruth_016"
save(niter, truemu, trueSigma, SXA, SXB, SIGMA_X0,
     KAPPA, TAU, Tx, Ty, idx, tx, idy, ty, X, M,
     Y, BETA, GAMMA, SIGMA_B, SIGMA_Y, ALPHA,
     file=file.path(paste0("samples/", filename,".RData")))
saveRDS(samples, file=file.path(paste0("samples/", filename,".RDS")))
print(paste0("Done, saved in file ", filename))


