Mufa <- function(niter, data, nthin=1, nburn=10, nchain=2,
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
  for (i in 1:nchain) {
    ### For Regression
    Y <- cst$Y
    prm[[i]][["B"]]        <- matrix(rnorm(cst$Kn*cst$q, 0, 10), cst$Kn, cst$q)
    prm[[i]][["Sigma"]]    <- diag(1/rgamma(cst$q, 2, 1), cst$q) 
    prm[[i]][["Sigmainv"]] <- solve(prm[[i]]$Sigma)
    prm[[i]][["l"]]        <- rep(rgamma(cst$q, 1, 0.1), cst$q)            # hierarchical IW
    prm[[i]][["nu"]]       <- rep(cst$aB*cst$tauB, cst$Kn)                 # shrinkage on B
    prm[[i]][["zeta"]]     <- cst$uB*prm[[i]]$nu                                # shrinkage on B
    prm[[i]][["gamma"]]    <- sample(0:1, size=nrow(prm[[i]]$B), replace=TRUE)  # select rows of B
    prm[[i]][["pi_gamma"]] <- runif(1, 0, 1)
    ### For Factor model
    X <- data$X
    nx <- nrow(X); px <- ncol(X)
    L <- cst$L; K <- cst$K
    prm[[i]][["sigmax"]] <- diag(1/rgamma(px, 2, 1), px)
    prm[[i]][["theta"]]  <- array(rnorm(px*L, 0, 10), dim = c(px, L))  # L 1st eigenvectors as initial conditions
    prm[[i]][["xi"]]     <- array(rnorm(nx * L * K), dim = c(nx, L, K))
    prm[[i]][["psi"]]    <- array(rnorm(nx * K), dim = c(nx, K))
    prm[[i]][["tau"]]    <- 1      # scale param for psi and xi, for identifiability, set to 1 while kappa varies freely
    prm[[i]][["eta"]]    <- array(rnorm(nx*K, 0, 10), dim = c(nx, K))
    prm[[i]][["phi"]]    <- array(rgamma(px * L, 1), dim = c(px, L))    # Shrinkage on theta
    prm[[i]][["delta"]]  <- rgamma(L, 1)                                    # Shrinkage on theta
    #----------------------
    # TODO: estimate from data bandwidth parameter for psi (and xi??)
    # model seems pretty robust to the choice of KAPPA
    if (is.null(data$KAPPA)) {
      prm[[i]][["kappa"]]  <- 10  
    } else {
      prm[[i]][["kappa"]]  <- data$KAPPA  
    }
    #----------------------
    prm[[i]][["Kgp"]]    <- se_cov(cst$tx, sigma = prm[[i]]$tau, kappa = prm[[i]]$kappa)  # TODO: estmate tau and kappa
    prm[[i]][["Kinv"]]   <- FastGP::rcppeigen_invert_matrix(prm[[i]]$Kgp)
    prm[[i]][["acp"]]    <- rep(0, nrow(prm[[i]]$eta))  # store acceptance rate for eta HM 
  }
  
  # Storage for outputs
  ##########################################################
  nout <- floor((niter-nburn)/nthin)
  out <- list()
  for (i in 1:nchain) {
    out[[i]] <- list(
      Lambda = array(NA, dim = c(nout, nx, px, K)),
      Mu = array(NA, dim = c(nout, nx, px)),
      Psi = array(NA, dim = c(nout, nx, px, px)),  # Cov(X)
      eta = array(NA, dim = c(nout, nx, K)),
      B     = array(NA, dim = c(nout, cst$Kn, cst$q)),
      Sigma = array(NA, dim = c(nout, cst$q, cst$q)), # Cov(Y)
      NULL) 
  }
  
  
  # Sampling
  #########################################################
  pb <- txtProgressBar(style = 3)
  j <- 1
  mh_delta <- rep(mh_delta, nchain)
  for (i in 1:niter) {
    for (c in 1:nchain) {
      ### Covariance regression ###
      prm[[c]] <- s_xi(prm[[c]], cst)  # time elapsed for 10^2 iters: 50.812 *
      prm[[c]] <- s_psi(prm[[c]], cst)  # time elapsed for 10^2 iters: 26.589 *
      prm[[c]] <- s_theta(prm[[c]], cst)  # time elapsed for 10^2 iters: 0.449 
      prm[[c]] <- s_sigmax(prm[[c]], cst)  # time elapsed for 10^2 iters: 5.878
      prm[[c]] <- s_eta(prm[[c]], cst, mh_delta[c])  # time elapsed for 10^2 iters: ?? 
      ### Y Regression ###
      if (method == "shrink") {
        stopifnot(sum(prm[[c]]$gamma)==cst$Kn)
        prm[[c]] <- s_B(prm[[c]], cst)
      } else if (method == "select") {
        prm[[c]] <- s_Bsns(prm[[c]], cst)
        
      } else stop(paste("Method", method, "not implemented"))
      prm[[c]] <- s_Sigma(prm[[c]], cst)
      
      # Update step_size for s_eta
      if (i%%50==0 & i<=nburn) {
        acp_mean = mean(prm[[c]]$acp)/100
        if(acp_mean > 0.45) {
          mh_delta[c] = mh_delta[c] + min(0.01, i^(-0.5))
        } else if(acp_mean < 0.2) {
          mh_delta[c] = mh_delta[c] - min(0.01, i^(-0.5))
        }
        prm[[c]][["acp"]] = rep(0, nrow(prm[[c]]$eta))
      }
      
      setTxtProgressBar(pb, i / niter)
      
      if (i>nburn & i%%nthin==0) {  # TODO: Is this a good way to thin????
        out[[c]]$eta[j,,] <- prm[[c]]$eta
        out[[c]]$Lambda[j,,,] <- abind::abind(lapply(1:nx, function(i) prm[[c]]$theta%*%prm[[c]]$xi[i,,]), along=-1)
        out[[c]]$Mu[j,,]      <- abind::abind(lapply(1:nx, function(i) out[[c]]$Lambda[j,i,,]%*%prm[[c]]$psi[i,]), along=-1) 
        out[[c]]$Psi[j,,,]    <- abind::abind(lapply(1:nx, function(i) out[[c]]$Lambda[j,i,,]%*%t(out[[c]]$Lambda[j,i,,]) + diag(prm[[c]]$sigmax)), along=-1)
        
        out[[c]]$B[j,,] <- prm[[c]]$B
        out[[c]]$Sigma[j,,] <- prm[[c]]$Sigma
        
        j <- j+1
      } 
    }
  }
  
  return(out)
}

