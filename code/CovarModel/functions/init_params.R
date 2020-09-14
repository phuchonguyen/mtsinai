# TODO: change sx to sxinv in all s_*.R files
# TODO: Let L or K varies ???
# TODO: add NA and LOD

init_params <- function(X, Y, 
                        ty, idy,
                        tx, idx, # for GP covariance of Y model
                        tau, kappa,  # TODO: estimate these
                        K, L  # TODO: Allow this to vary
                        ) {
  # X = rbind(X_1,...,X_n), X_1 = rbind(X_11,...X_1Tx)
  temp <- data.frame(cbind(X, idx, tx)) %>%
    arrange(idx, tx) 
  X <- temp %>%
    select(-idx, -tx) %>% as.matrix()
  idx <- temp %>% select(idx) %>% as.matrix() %>% as.vector()
  tx <- temp %>% select(tx) %>% as.matrix() %>% as.vector()
  
  Tx <- length(unique(tx))
  
  # X latent variable initial values
  nx <- nrow(X); px <- ncol(X)
  eta <- array(rnorm(nx*K), dim = c(nx, K))
  xi <- array(rnorm(nx*L*K), dim = c(nx, L, K))
  psi <- array(rnorm(nx*K), dim = c(nx, K))
  theta <- array(rnorm(px*L), dim = c(px, L))
  sigmax <- 1/rgamma(px, 1, 1)
  phi <- array(rgamma(px*L, 1), dim = c(px, L))
  delta <- rgamma(L, 1)
  Kgp <- se_cov(tx, sigma = tau, kappa = kappa)  # TODO: estmate tau and kappa
  Kinv <- rcppeigen_invert_matrix(Kgp) # TODO: optimize?
  
  params <- list(
    #----------------------------
    # For covariance factor model
    #---------------------------
    theta = theta,
    xi = xi,
    sigmax = sigmax,
    psi = psi,
    kappa = kappa,  # bandwidth parameter for psi (and xi??)
    tau = tau,      # scale param for psi and xi ???
    eta = eta,
    Kinv = Kinv,    # TODO: update if kappa and tau are updated.
    phi = phi,
    delta = delta
  )
  
  const <- list(
    # For covariance factor model
    X = X,
    L = L,  # TODO: how to not fix these?
    K = K,  # TODO: let K varies?
    tx = tx,
    idx = idx,
    Tx = Tx, 
    a1delta = 2,  # params for shrinkage factor, from simulation in Fox 2015
    a2delta = 2,  # params for shrinkage factor, from simulation in Fox 2015
    bphi = 3,     # params for shrinkage factor, from simulation in Fox 2015
    asx = 2,      # for Gamma prior on Sigma_X = diag(sigmax_1,...,sigmax_p), value from simulation
    bsx = 0.1,    # for Gamma prior on Sigma_X, value from simulation
    asy = 1,      # for Gamma prior on Sigma_Yt
    bsy = 1,      # for Gamma prior on Sigma_Yt
    asb = 1,      # for Gamma prior on Sigma_beta
    bsb = 1      # for ''  Sigma_beta
  )
  
  return(list(params=params, const=const))
}