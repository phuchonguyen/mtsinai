# TODO clean up

s_xi <- function(prm, cst, args=NULL) {
  L <- cst$L
  K <- cst$K
  X <- cst$X
  theta <- prm$theta
  xi <- prm$xi
  eta <- prm$eta
  sigmax <- prm$sigmax
  Kinv <- prm$Kinv
  n <- nrow(X)
  
  #sxinv <- diag(1/prm$sigmax) # Sigma_0 is a diagonal matrix
  #S <- bdiag(lapply(1:n, function(i) sxinv))  # Sparse matrix 
  mu <- lapply(1:n, function(i) X[i,] - theta%*%xi[i,,]%*%eta[i,])
  
  for (l in 1:L) {
    for (k in 1:K) {
      Xtilda <- lapply(1:n, function(i) mu[[i]] + eta[i,k]*xi[i,l,k]*theta[,l])
      #Xtilda <- unlist(Xtilda)  # stack column vectors Xi
      #A <- get_A(X, eta, xi, theta, l, k)
      #AtS <- t(A)%*%S
      #Stilda <- rcppeigen_invert_matrix(Kinv + as.matrix(AtS%*%A)) #TODO: Optimize inversion here, convert sparse to dense  
      #mutilda <- as.matrix(Stilda%*%AtS%*%Xtilda)  # convert sparse to dense
      sn <- sapply(1:n, function(i) eta[i,k]^2*sum(theta[,l]^2/sigmax))
      Stilda <- rcppeigen_invert_matrix(Kinv + diag(sn))
      mutilda <- sapply(1:n, function(i) eta[i,k]*sum(theta[,l]*Xtilda[[i]]/sigmax))
      xi[,l,k] <- rcpp_rmvnorm(n = 1, S = Stilda, mu = Stilda%*%mutilda)
    }
  }
  
  prm[["xi"]] <- xi
  return(prm)
}


# Returns a sparse matrix
get_A <- function(X, eta, xi, theta, l, k) {
  # n <- nrow(X)
  # p <- ncol(X)
  # val <- as.vector(t(eta[,k]%*%t(theta[,l])))
  # indc <- rep(1:n, each=p)
  # indr <- 1:(n*p)
  # A <- sparseMatrix(
  #   i = indr, 
  #   j = indc, 
  #   x = val,
  #   dims = c(n*p, n))
  eta_k <- eta[,k]
  theta_l <- theta[,l]
  A <- diag(eta_k) %x% theta_l
  return(A)
}
