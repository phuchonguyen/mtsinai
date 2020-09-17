# TODO: change sx to sxinv in all s_*.R files
# TODO: Let L or K varies ???
# TODO: add NA and LOD

init_params <- function(X, Y, 
                        ty, idy,
                        K  # TODO: Allow this to vary
                        ) {
  
  # Y = c(Y_1,...,Y_Ty), Y_1 = c(y_11,...,y_1n)
  temp <- data.frame(cbind(Y, idy, ty)) %>%
    arrange(ty, idy) 
  Y <- temp %>%
    select(Y) %>% as.matrix() %>% as.vector()
  idy <- temp %>% select(idy) %>% as.matrix() %>% as.vector()
  ty <- temp %>% select(ty) %>% as.matrix() %>% as.vector()
  Ty <- length(unique(ty))
  
  # Y regression initial values
  py <- ncol(X)
  pi_gamma <- 0.5
  gamma <- rbinom(py, 1, pi_gamma)
  beta <- matrix(0, Ty, py)
  alpha <- sapply(1:Ty, function(t) mean(Y[ty==t]))
  #beta_int <- rep(0,p)
  #beta_z <- rep(0,ncol(Z))  #TODO: add this
  sigmay <- rep(1, Ty)
  sigmab <- rep(1, py)
  V <- rw_cov(Ty)  # TODO generalize to other model for beta_t??? This V is unimodular!
  Vdet <- det(V)
  Vinv <- solve(V)
  
  
  params <- list(
    #-----------------------
    # Y regression
    #----------------------
    pi_gamma = pi_gamma,
    gamma = gamma,
    beta = beta,
    alpha = alpha,
    #beta_int = beta_int,
    #beta_z = beta_z,  # TODO: Add
    sigmay = sigmay,
    sigmab = sigmab
  )
  
  const <- list(
    # For covariance factor model
    X = X,
    K = K,  # TODO: let K varies?
    Y = Y,
    ty = ty,
    idy = idy,
    Ty = Ty,
    asb = 1,
    bsb = 0.1,
    asy = 1,      # for Gamma prior on Sigma_Yt
    bsy = 1,      # for Gamma prior on Sigma_Yt
    V = V,        # Covariance of beta_.j
    Vdet = Vdet,
    Vinv = Vinv,
    pi_gamma0 = 0.5
  )
  
  return(list(params=params, const=const))
}