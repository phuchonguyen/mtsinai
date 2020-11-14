Mufa <- function(niter, data, nthin=1, nburn=10,
                 mh_delta=0.05, method="shrink") {
  
  # Constants
  cst <- list(X = data$X,        # (total # measurements x # chemicals)
              tx = data$tx,      # time of each exposure measurement
              idx = data$idx,    # subject ID's of each exposure measurement
              Tx = max(data$tx),
              p = ncol(data$X),
              K = data$K,
              L = data$L,
              Y = data$Y,        # (# subjects x # outcomes), assumes Y is sorted by increasing subject ID
              q = ncol(data$Y),  # m in model
              N = nrow(data$Y),  # number of subjects
              Kn = data$K*max(data$tx),
              NULL)

  # Setting hyper params 
  ############################################################
  ### Hyperparams for Factor Model
  cst[["a1delta"]] <- 2  # factor global shrinkage 
  cst[["a2delta"]] <- 2  # factor global shrinkage
  if (is.null(data$aphi)) {
    cst[["aphi"]]    <- matrix(1, cst$p, cst$L)     # factor local shrinkage
    cst[["bphi"]]    <- matrix(3, cst$p, cst$L)     # factor local shrinkage
  } else {
    stopifnot(dim(data$aphi) == c(cst$p, cst$L))
    stopifnot(dim(data$bphi) == c(cst$p, cst$L))
    cst[["aphi"]]    <- data$aphi     # factor local shrinkage
    cst[["bphi"]]    <- data$bphi     # factor local shrinkage
  }
  cst[["asx"]]     <- 2      # for Gamma prior on Sigma_X = diag(sigmax_1,...,sigmax_p)
  cst[["bsx"]]     <- 0.1    # for Gamma prior on Sigma_X
  ### Hyperparams for Regression Model
  cst[["vS"]] <- 2  # param in prior of Sigma
  cst[["dS"]] <- rep(10, cst$q)  # hyperparam in prior of Sigma
  cst[["uB"]] <- 1/2  # hyper param for shrinkage on B
  cst[["aB"]] <- 1/2
  cst[["tauB"]] <- 1/(cst$q*sqrt(cst$N * log(cst$N)))  # as recommended by Bai & Ghosh
  cst[["api"]] <- 1   # hyperparam for pi_gamma
  cst[["bpi"]] <- 1   # hyperparam for pi_gamma
  

  # Initial conditions
  #########################################################
  prm <- list()
  ### For Regression
  Y <- cst$Y
  prm[["B"]]        <- matrix(1, cst$Kn, cst$q)
  prm[["Sigma"]]    <- diag(1, cst$q) 
  prm[["Sigmainv"]] <- solve(prm$Sigma)
  prm[["l"]]        <- rep(1, cst$q)                 # hierarchical IW
  prm[["nu"]]       <- rep(cst$aB*cst$tauB, cst$Kn)  # shrinkage on B
  prm[["zeta"]]     <- cst$uB*prm$nu                 # shrinkage on B
  prm[["gamma"]]    <- rep(1, nrow(prm$B))           # select rows of B
  prm[["pi_gamma"]] <- 0.5
  ### For Factor model
  X <- data$X
  nx <- nrow(X); px <- ncol(X)
  L <- cst$L; K <- cst$K
  prm[["sigmax"]] <- apply(X, 2, var)
  svd_out <- svd(X - colMeans(X))
  prm[["theta"]]  <- svd_out$v[,1:L]  # L 1st eigenvectors as initial conditions
  prm[["xi"]]     <- array(rnorm(nx * L * K), dim = c(nx, L, K))
  prm[["psi"]]    <- array(rnorm(nx * K), dim = c(nx, K))
  prm[["tau"]]    <- 1      # scale param for psi and xi, for identifiability, set to 1 while kappa varies freely
  prm[["eta"]]    <- svd_out$u[,1:K] %*% diag(svd_out$d[1:K])
  prm[["phi"]]    <- array(rgamma(px * L, 1), dim = c(px, L))    # Shrinkage on theta
  prm[["delta"]]  <- rgamma(L, 1)                                    # Shrinkage on theta
  #----------------------
  # TODO: estimate from data bandwidth parameter for psi (and xi??)
  # model seems pretty robust to the choice of KAPPA
  if (is.null(data$KAPPA)) {
    prm[["kappa"]]  <- 10  
  } else {
    prm[["kappa"]]  <- data$KAPPA  
  }
  #----------------------
  prm[["Kgp"]]    <- se_cov(cst$tx, sigma = prm$tau, kappa = prm$kappa)  # TODO: estmate tau and kappa
  prm[["Kinv"]]   <- FastGP::rcppeigen_invert_matrix(prm$Kgp)
  prm[["acp"]]    <- rep(0, nrow(prm$eta))  # store acceptance rate for eta HM
  
  # Storage for outputs
  ##########################################################
  nout <- floor((niter-nburn)/nthin)
  out <- list(
    Lambda = array(NA, dim = c(nout, nx, px, K)),
    Mu = array(NA, dim = c(nout, nx, px)),
    Psi = array(NA, dim = c(nout, nx, px, px)),  # Cov(X)
    eta = array(NA, dim = c(nout, nx, K)),
    B     = array(NA, dim = c(nout, cst$Kn, cst$q)),
    Sigma = array(NA, dim = c(nout, cst$q, cst$q)), # Cov(Y)
    NULL)
  
  
  # Sampling
  #########################################################
  pb <- txtProgressBar(style = 3)
  j <- 1
  for (i in 1:niter) {
    ### Covariance regression ###
    prm <- s_xi(prm, cst)  # time elapsed for 10^2 iters: 50.812 *
    prm <- s_psi(prm, cst)  # time elapsed for 10^2 iters: 26.589 *
    prm <- s_theta(prm, cst)  # time elapsed for 10^2 iters: 0.449 
    prm <- s_sigmax(prm, cst)  # time elapsed for 10^2 iters: 5.878
    prm <- s_eta(prm, cst, mh_delta)  # time elapsed for 10^2 iters: ?? 
    ### Y Regression ###
    if (method == "shrink") {
      stopifnot(sum(prm$gamma)==cst$Kn)
      prm <- s_B(prm, cst)
    } else if (method == "select") {
      prm <- s_Bsns(prm, cst)
      
    } else stop(paste("Method", method, "not implemented"))
    prm <- s_Sigma(prm, cst)
    
    # Update step_size for s_eta
    if (i%%50==0 & i<=nburn) {
      acp_mean = mean(prm$acp)/100
      if(acp_mean > 0.45) {
        mh_delta = mh_delta + min(0.01, i^(-0.5))
      } else if(acp_mean < 0.2) {
        mh_delta = mh_delta - min(0.01, i^(-0.5))
      }
      prm[["acp"]] = rep(0, nrow(prm$eta))
    }
    
    setTxtProgressBar(pb, i / niter)
    
    if (i>nburn & i%%nthin==0) {  # TODO: Is this a good way to thin????
      out$eta[j,,] <- prm$eta
      out$Lambda[j,,,] <- abind::abind(lapply(1:nx, function(i) prm$theta%*%prm$xi[i,,]), along=-1)
      out$Mu[j,,]      <- abind::abind(lapply(1:nx, function(i) out$Lambda[j,i,,]%*%prm$psi[i,]), along=-1) 
      out$Psi[j,,,]    <- abind::abind(lapply(1:nx, function(i) out$Lambda[j,i,,]%*%t(out$Lambda[j,i,,]) + diag(prm$sigmax)), along=-1)
      
      out$B[j,,] <- prm$B
      out$Sigma[j,,] <- prm$Sigma
      
      j <- j+1
    }
  }
  
  return(out)
}

