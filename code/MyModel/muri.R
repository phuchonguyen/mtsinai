Muri <- function(niter, data, method="shrink") {
  cst <- list(Y = data$Y,  # N x m matrix in model
              X = data$X,
              Tx = 1,             # TODO Update # n in model  
              q = ncol(data$Y),  # m in model
              N = nrow(data$Y),  # number of subjects
              Kn = ncol(data$X),
              NULL)
  
  # Setting hyperparam
  cst[["vS"]] <- 2  # param in prior of Sigma
  cst[["dS"]] <- rep(10, cst$q)  # hyperparam in prior of Sigma
  cst[["uB"]] <- 1/2  # hyper param for shrinkage on B
  cst[["aB"]] <- 1/2
  cst[["tauB"]] <- 1/(cst$q*sqrt(cst$N * log(cst$N)))  # as recommended by Bai & Ghosh
  cst[["api"]] <- 1   # hyperparam for pi_gamma
  cst[["bpi"]] <- 1   # hyperparam for pi_gamma
  
  # Setting initial conditions
  X <- cst$X
  Y <- cst$Y
  XtX <- t(X)%*%X
  XtY <- t(X)%*%Y
  min.sing <- min(svd(XtX)$d)
  delta <- 0.01
  prm <- list()   
  prm[["B"]]        <- chol2inv(chol(XtX + (delta+min.sing)*diag(cst$Kn)))%*% XtY# Initial guess MLE from MBSP package
  prm[["Sigma"]]    <- (n-1)/n * cov(Y - X%*%prm$B)  # Initial guess MLE from MBSP package
  prm[["Sigmainv"]] <- solve(prm$Sigma)
  prm[["l"]]        <- rep(1, cst$q)  # hierarchical IW
  prm[["nu"]]       <- rep(cst$aB*cst$tauB, cst$Kn)  # shrinkage on B
  prm[["zeta"]]     <- cst$uB*prm$nu                 # shrinkage on B
  prm[["gamma"]]    <- rep(1, nrow(prm$B))           # select rows of B
  prm[["pi_gamma"]] <- 0.5
  
  # Storage for outputs
  out <- list(B     = array(NA, dim = c(niter, cst$Kn, cst$q)),
              Sigma = array(NA, dim = c(niter, cst$q, cst$q)),
              NULL)
  
  pb = txtProgressBar(style = 3)
  for (i in 1:niter) {
    
    ### Sample B ###
    if (method == "shrink") {
      stopifnot(sum(prm$gamma)==cst$Kn)
      prm <- s_B(prm, cst)
      
    } else if (method == "select") {
      prm <- s_Bsns(prm, cst)
    
    } else stop(paste("Method", method, "not implemented"))
    
    ### Sample Sigma ###
    prm <- s_Sigma(prm, cst)
    
    setTxtProgressBar(pb, i / niter)
    
    out$B[i,,] <- prm$B
    out$Sigma[i,,] <- prm$Sigma
  }
  
  return(out)
}