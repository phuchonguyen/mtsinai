# TODO: change sx to sxinv in all s_*.R files
# TODO: Let L or K varies ???
# TODO: add NA and LOD

init_params <- function(X, Y, 
                        ty, idy,
                        tx, idx, # for GP covariance of Y model
                        tau, kappa,  # TODO: estimate these
                        K, L) {
  
  # Scale time to be in interval [0, 1]
  if (scale.t) {
    ty <- (ty - min(ty))/(max(ty) - min(ty))
    tx <- (tx - min(tx))/(max(tx) - min(tx))
  }
  
  # Y = c(Y_1,...,Y_Ty), Y_1 = c(y_11,...,y_1n)
  temp <- data.frame(cbind(Y, idy, ty)) %>%
    arrange(ty, idy) 
  Y <- temp %>%
    select(Y) %>% as.matrix() %>% as.vector()
  idy <- temp %>% select(idy) %>% as.matrix() %>% as.vector()
  ty <- temp %>% select(ty) %>% as.matrix() %>% as.vector()
  
  # X = rbind(X_1,...,X_n), X_1 = rbind(X_11,...X_1Tx)
  temp <- data.frame(cbind(X, idx, tx)) %>%
    arrange(idx, tx) 
  X <- temp %>%
    select(-idx, -tx) %>% as.matrix()
  idx <- temp %>% select(idx) %>% as.matrix() %>% as.vector()
  tx <- temp %>% select(tx) %>% as.matrix() %>% as.vector()
  
  Tx <- length(unique(tx))
  Ty <- length(unique(ty))
  
  # Y regression initial values
  py <- K*Tx
  pi_gamma <- 0.5
  gamma <- rbinom(py, 1, pi_gamma)
  beta <- matrix(0, Ty, py)
  beta0 <- rep(0, py)
  alpha <- sapply(1:Ty, function(t) mean(Y[ty==t]))
  #beta_int <- rep(0,p)
  #beta_z <- rep(0,ncol(Z))  #TODO: add this
  sigmay <- rep(1, Ty)
  sigmab <- rep(1, py)
  V <- rw_cov(Ty)  # TODO generalize to other model for beta_t??? This V is unimodular!
  Vdet <- det(V)
  Vinv <- solve(V)
  
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
    #-----------------------
    # Y regression
    #----------------------
    pi_gamma = pi_gamma,
    gamma = gamma,
    beta = beta,
    beta0 = beta0,
    alpha = alpha,
    #beta_int = beta_int,
    #beta_z = beta_z,  # TODO: Add
    sigmay = sigmay,
    sigmab = sigmab,
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
    acp = rep(0, nrow(eta)),
    eta_naccept = 0,
    eta_npro = 0,
    Kinv = Kinv,    # TODO: update if kappa and tau are updated.
    phi = phi,
    delta = delta
  )
  
  const <- list(
    # For covariance factor model
    X = X,
    L = L,  # TODO: how to not fix these?
    K = K,  # TODO: let K varies?
    Y = Y,
    ty = ty,
    tx = tx,
    idx = idx,
    idy = idy,
    Tx = Tx, 
    Ty = Ty,
    a1delta = 2,  # params for shrinkage factor, from simulation in Fox 2015
    a2delta = 2,  # params for shrinkage factor, from simulation in Fox 2015
    bphi = 3,     # params for shrinkage factor, from simulation in Fox 2015
    asx = 2,      # for Gamma prior on Sigma_X = diag(sigmax_1,...,sigmax_p), value from simulation
    bsx = 0.1,    # for Gamma prior on Sigma_X, value from simulation
    asy = 1,      # for Gamma prior on Sigma_Yt
    bsy = 1,      # for Gamma prior on Sigma_Yt
    V = V,        # Covariance of beta_.j
    Vdet = Vdet,
    Vinv = Vinv,
    asb = 1,      # for Gamma prior on Sigma_beta
    bsb = 1,      # for ''  Sigma_beta
    mb0 = rep(0, py),      # mean beta0
    sb0 = rep(10, py),      # variance beta0
    pi_gamma0 = 0.5   # prior prob. of gammaj = 0
  )
  
  return(list(params=params, const=const))
}