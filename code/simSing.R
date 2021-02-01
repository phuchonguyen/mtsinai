#.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))
source(file.path(getwd(), "code/MyModel/mufa.R"))
source(file.path(getwd(), "code/sim_func.R"))
funcpath <- file.path(getwd(), "code/MyModel/functions")
funcfiles <- list.files(path = funcpath)
sapply(funcfiles, function(x) source(file.path(funcpath, x)))

# Constants
SEED <- 123
set.seed(SEED)

# Model parameters
K <- 2
Tx <- 3
Ty <- 3
P <- 10
L <- 3
M <- 100    # Number of subjects
KAPPA <- 10 # GP bandwidth
EXP_P <- 2  # GP Gaussian kernel
TAU <- 2   # GP sqrt variance
SIGMA_X0 <- diag(0.5, P, P)#diag(1/rgamma(P, SXA, SXB), P, P)

# Generate X
# Generate Psi
truepsi <- array(NA, dim = c(Tx, K))
for (i in 1:K) {
  truepsi[,i] <- rgp(1, 1:Tx, mu_zero, KAPPA, TAU)
}
tx <- rep(1:Tx, M)  # TODO: Allow exposures to be missing at different times
psi <- t(sapply(tx, function(t) truepsi[t,]))
psi <- as.matrix(psi)
if (K == 1 & ncol(psi) != 1) {psi <- t(psi)}
eta <- t(apply(psi, 1, function(x) mvtnorm::rmvnorm(1, mean=x)))
eta <- as.matrix(eta)
if (K == 1 & ncol(eta) != 1) {eta <- t(eta)}
idx <- rep(1:M, each=Tx)
etay <- transform_etay(eta, idx, Tx)

# Generate Theta
a1 <- a2 <- gtheta <- 10
truedelta <- rgamma(1, a1, 1)
truedelta <- c(truedelta, rgamma(L-1, a2, 1))
truephi <- matrix(rgamma(L*P, gtheta/2, gtheta/2), P, L)
truetau <- get_factor_tau(truedelta)
Theta <- array(NA, dim = c(P, L))
for (j in 1:P) {
  for (l in 1:L) {
    Theta[j,l] <- rnorm(1, 0, 1/sqrt(truephi[j,l] * truetau[l]))
  }
}

# Generate Xi
truexi <- array(NA, dim = c(Tx, L, K))
for (k in 1:K) {
  for (l in 1:L) {
    truexi[,l,k] <- rgp(1, 1:Tx, mu_zero, KAPPA, TAU/2)
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
# Visualize Mu & Sigma
#my_plot_gp(mu = truemu, Sigma = trueSigma, Tx = Tx)

# Generate X
X <- t(sapply(1:nrow(eta), function(i) 
  mvtnorm::rmvnorm(1, 
                   mean=Theta%*%truexi[tx[i],,]%*%eta[i,],
                   sigma=SIGMA_X0)))

# Generate Y
q <- 3   # number of outcomes
n <- M
p <- K*Tx
# Active number of predictors is :
p.act <- 2
# Generate true error covariance Sigma Y
rho <- 0.6
sigma.sq <- 0.5
times <- 1:q
H <- abs(outer(times, times, "-"))
Sigma <- sigma.sq * rho^H
# Generate noise matrix E 
E <- MASS::mvrnorm(n, rep(0,q), Sigma)
B.act <- MASS::mvrnorm(p.act, seq(1, 2, length.out = q), Sigma)  # Correlated B
B.true <- rbind(B.act,matrix(0,p-p.act,q))
# Generate Y
Y <- etay%*%B.true + E


##########################################################
# Run Model  
##########################################################

# Running it for longer improves the estimates.
data <- list(X=X, tx=tx, idx=idx, K=K, L=L, KAPPA=KAPPA,
             #aphi=aphi, bphi=bphi,
             Y=Y)
s <- Sys.time()
samps <- Mufa(100, data, nburn = 20, nthin=2)
saveRDS(samps, paste("code/samples/sim5/sim005_K",K, "Tx", Tx, "Ty",q, Sys.Date(), ".RDS", sep="_"))
save(SIGMA_X0, truemu, trueSigma, eta,
     B.true, Sigma, Theta, truepsi, truexi, X, Y,
     file=paste("code/samples/sim5/sim005_K",K, "Tx", Tx, "Ty",q, Sys.Date(), ".Rdata", sep="_"))
print(Sys.time()-s)









