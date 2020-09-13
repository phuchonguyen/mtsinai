# Convert eta_i from "long" form from Factor Model to "wide" form for regression
# Ex: rbind(eta_i1, eta_i2, eta_i3) to c(eta_i1, eta_i2, eta_i3)
to_etai_reg <- function(etai_fac) {
  return(as.vector(t(etai_fac)))
}


# Generate covariance matrix for Squared Exponential kernel
# c(x, x') = sigma^2 exp(-kappa||x-x'||^2_2)4
# Add a perturbation to make positive definite
# x: array of points (n x p)
# Returns:
# (n x n) positive-definite matrix
se_cov <- function(x, sigma, kappa) {
  # If p = 1
  p <- ncol(as.matrix(x))
  if (p==1) {
    x <- as.vector(x)
    perturb <- diag(rep(1e-10, length(x)))
    return(outer(x, x, function(a, b) sigma^2*exp(-1/kappa*(a-b)^2)) + perturb)
  }
  perturb <- diag(rep(1e-10, nrow(x)))
  if (length(kappa) < p) kappa <- rep(kappa, ncol(x))
  x <- sweep(x, 2, sqrt(1/kappa), "*")
  K <- sigma^2*fields::Exp.cov(x, theta=1, p=2) # p=2 is Squared Exponential Function
  return(K + perturb)
}


# TESTED!
# Make random walk covariance matrix
rw_cov <- function(p) {
  m <- matrix(NA, p, p)
  for (j in 1:p) {
    m[j:p, j:p] <- j
  }
  return(m)
}


# Calculate mean vector for y
# j : remove the jth covariate from the mean
# TODO add quadratic term
get_muy <- function(params, const, eta=NULL, j=NULL) {
  if (is.null(eta)) {
    eta <- params$eta
  }
  etay <- transform_etay(eta, const$idx, const$Tx)
  #etay_int <- etay^2
  if (is.null(j)) {
    muy <- sapply(1:length(Y), function(i) {
      const$Y[i] - t(etay[const$idy[i],])%*%(params$beta[const$ty[i], ]) - params$alpha[const$ty[i]]
    }) 
  } else {
    muy <- sapply(1:length(Y), function(i) {
      const$Y[i] - t(etay[const$idy[i],-j])%*%(params$beta[const$ty[i], -j]) - params$alpha[const$ty[i]]
    })
  }
  #+etay_int[const$idy,]%*%params$beta_int)
  return(muy)
}


# Returns: a (n x KTx) matrix with vectorized exposures of all time for each individual
# So each individual has one row, sorted by increasing ID
transform_etay <- function(eta, idx, Tx) {
  uidx <- sort(unique(idx))
  etay <- t(sapply(uidx, function(i) as.vector(t(eta[idx==i,]))))
  eta_names <- sapply(1:ncol(eta), function(k) paste("eta", k, sep=""))
  etay_names <- as.vector(t(sapply(1:Tx, function(t) paste(eta_names, t, sep="_"))))
  colnames(etay) <- etay_names
  rownames(etay) <- uidx
  return(etay)
}


# Orthogonal projection onto col(X)
get_projection_matrix <- function(X) {
  n <- nrow(X)
  if (is.null(n)) {n <- 1}
  return(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))
}
