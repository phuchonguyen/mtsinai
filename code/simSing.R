.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))
source(file.path(getwd(), "MyModel/mufa.R"))
source(file.path(getwd(), "sim_func.R"))
funcpath <- file.path(getwd(), "MyModel/functions")
funcfiles <- list.files(path = funcpath)
sapply(funcfiles, function(x) source(file.path(funcpath, x)))

# Constants
SEED <- 123
set.seed(SEED)

K <- 2 
Tx <- 3
P <- 10
L <- 3
M <- 50    # Number of subjects
KAPPA <- 100 # GP bandwidth
EXP_P <- 2  # GP Gaussian kernel
TAU <- 1    # GP sqrt variance
SXA <- 2
SXB <- 0.1
SIGMA_X0 <- diag(1/rgamma(P, SXA, SXB), P, P)  # Error variance of X
tx <- rep(1:Tx, M)  # TODO: Allow exposures to be missing at different times
idx <- rep(1:M, each=Tx)

# Generate Theta
# -----------Option 1 -----------
a1 <- a2 <- gtheta <- 10
truedelta <- rgamma(1, a1, 1)
truedelta <- c(truedelta, rgamma(L-1, a2, 1))
truephi <- matrix(rgamma(L*P, gtheta/2, gtheta/2), P, L)
# ------------Option 2: more sparse Theta----------------
#truedelta <- c(rgamma(1, 2, 1), rgamma(1, 3, 1), rgamma(L-2, 10, 1))
#aphi <- bphi <- matrix(3, P, L)
#for (l in 1:L) {
#  s <- 1+(l-1)*ceiling(P/L)
#  e <- min(P, l*ceiling(P/L))
#  aphi[s:e, l] <- 1
#  aphi[-(s:e), l] <- 6
#}
#truephi <- t(matrix(sapply(1:P, function(i) rgamma(L, aphi[i]/2, bphi[i]/2)), L, P))
truetau <- get_factor_tau(truedelta)
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
# Calculate Mu & Sigma
truemu <- array(NA, dim = c(Tx, P))
trueSigma <- array(NA, dim = c(Tx, P, P))
for (t in 1:Tx) {
  Lambda <- Theta%*%truexi[t,,]
  truemu[t,] <- Lambda%*%truepsi[t,]
  trueSigma[t,,] <- Lambda%*%t(Lambda) + SIGMA_X0
}
#' ### Generate X from the truth
# Generate eta
psi <- t(sapply(tx, function(t) truepsi[t,]))
psi <- as.matrix(psi)
if (K == 1 & ncol(psi) != 1) {psi <- t(psi)}
eta <- t(apply(psi, 1, function(x) mvtnorm::rmvnorm(1, mean=x)))
eta <- as.matrix(eta)
if (K == 1 & ncol(eta) != 1) {eta <- t(eta)}
X <- t(sapply(1:nrow(eta), function(i) mvtnorm::rmvnorm(1, mean=Theta%*%truexi[tx[i],,]%*%eta[i,],
                                                        sigma=SIGMA_X0)))

##########################################################
# Generate Y
##########################################################
q <- 3   # number of outcomes
n <- M
p <- K*Tx
# Active number of predictors is :
p.act <- 1
# Generate true error covariance Sigma Y
rho <- 0.6
sigma.sq <- 2
times <- 1:q
H <- abs(outer(times, times, "-"))
Sigma <- sigma.sq * rho^H
# Generate noise matrix E 
mu <- rep(0,q)
E <- MASS::mvrnorm(n, mu, Sigma)
# Generate true coefficient matrix B_0. #
# Entries in nonzero rows are drawn from Unif[(-5,-0.5)U(0.5,5)]
# Option 1: B[i,j] iid from U(-5, 4)
#B.act <- runif(p.act*q,-5,4)
# Option 2: B[i,] iid from N(m, Sigma), E[i,] ~ N(0, Sigma)
B.act <- MASS::mvrnorm(p.act, seq(-5, 4, length.out = q), Sigma)  # Correlated B
disjoint <- function(x){
  if(x <= -0.5) return(x)
  else return(x+1)
}
B.act <- matrix(sapply(B.act, disjoint),p.act,q)
# Set rest of the rows equal to 0
B.true <- rbind(B.act,matrix(0,p-p.act,q))
B.true <- B.true[sample(1:p),] # permute the rows

# Generate response matrix Y #
etay <- transform_etay(eta, idx, Tx)
Y <- crossprod(t(etay),B.true) + E


##########################################################
# Run Model  
##########################################################

# Running it for longer improves the estimates.
data <- list(X=X, tx=tx, idx=idx, K=K, L=L, KAPPA=KAPPA,
             #aphi=aphi, bphi=bphi,
             Y=Y)
s <- Sys.time()
samps <- Mufa(60000, data, nburn = 20000, nthin=40)
saveRDS(samps, paste("samples/remote_K",K, "Tx", Tx, "Ty",q, Sys.Date(), ".RDS", sep="_"))
save(SIGMA_X0, truemu, trueSigma, eta,
     B.true, Sigma, Theta, truepsi, truexi, X, Y,
     file=paste("samples/remote_K",K, "Tx", Tx, "Ty",q, Sys.Date(), ".Rdata", sep="_"))
print(Sys.time()-s)









