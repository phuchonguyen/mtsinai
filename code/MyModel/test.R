X <- matrix(rnorm(100), 25, 4)
Px <- X%*%solve(t(X)%*%X)%*%t(X)
#Xsingular <- cbind(X, X)
#Pxsingular <- Xsingular%*%matlib::Ginv(t(Xsingular)%*%Xsingular)%*%t(Xsingular)
L <- 2; K <- 2
theta <- matrix(rnorm(4*L), 4, L)
xi <- array(1, dim=c(25,L,K))
eta <- matrix(1, 25, K)
psi <- matrix(1, 25, K)
sx <- diag(1, 4, 4)
delta <- 1:L
phi <- matrix(rnorm(4*L), 4, L)
kappa <- 2
a0 <- 1
b0 <- 1
Tx <- 5
Ty <- 3
ty <- 
tx <- rep(1)
Y <- const$Y
X <- const$X
idx <- const$idx
idy <- const$idy
eta <- params$eta
beta <- params$beta
beta_int <- params$beta_int
CI_inv <- params$CI_inv
theta <- params$theta
xi <- params$xi
psi <- params$psi
sigmax <- params$sigmax
