# Code from MixSelect

# Spike and slab prior for all non-linear terms
# This is the tau param in paper

s_rho = function(params, const) {
  
  rho = params$rho
  sigmay = params$sigmay
  
  
  # TODO: Check these 2
  X = transform_etay(params$eta, const$idx, const$Tx)
  mu = get_muy(params, const)
  p = ncol(X)
  n = nrow(X)
  
  
  P = params$P; Pt = t(P)
  PCPt = params$PCPt
  logdetCI = params$logdetCI
  CI_inv = params$CI_inv
  CI = params$CI
  exp_C = const$exp_C
  
  if (rho == 0){
    
    rho_star = rgamma(1,1)
    jstar = sample(1:p,1)
    lj_star = rgamma(1,1)
    l_star = rep(0,p)
    l_star[jstar] = lj_star
    X_sweep_star = sweep(X, 2, l_star, `*`)
    
    # covariance function
    C_star = rho_star*fields::Exp.cov(X_sweep_star, x2=NULL, p = exp_C)
    PCPt_star = P%*%C_star%*%Pt
    CI_star = PCPt_star + sigmay*diag(n)
    chol = chol(CI_star)
    CI_inv_star = chol2inv(chol)
    logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
    
    # TODO CHECK THIS PLZZZ
    loglik = logdetCI_star - logdetCI +
      sum(sapply(1:Ty, function(i) 
        t(mu[ty==i])%*%(CI_inv_star - CI_inv)%*%mu[ty==i]))
    loglik = -0.5*loglik
    
    logr = loglik
    logu = log(runif(1))
    
    if (logr > logu){
      rho = rho_star
      CI = (t(CI_star) + CI_star)/2
      CI_inv = (CI_inv_star + t(CI_inv_star))/2
      logdetCI = logdetCI_star
      PCPt = PCPt_star
      
      l = l_star
      X_sweep = X_sweep_star
      gamma_l = rep(0,p)
      gamma_l[jstar] = 1
      
    }else{
      l = rep(0,p)
      gamma_l = rep(0,p)
      X_sweep = sweep(X, 2, l, `*`)
    }
  }
  else{
    
    rho_star = 0
    
    PCPt_star = matrix(0,n,n)
    CI_star = diag(n)*sigmay
    chol = chol(CI_star)
    logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
    CI_inv_star = chol2inv(chol)
    
    # TODO: CHECK THIS!!!
    loglik = logdetCI_star - logdetCI +
      sum(sapply(1:Ty, function(i) 
        t(mu[ty==i])%*%(CI_inv_star - CI_inv)%*%mu[ty==i]))
    loglik = -0.5*loglik
    
    logr = loglik
    logu = log(runif(1))
    
    if (logr > logu){
      
      rho = rho_star
      CI = CI_star
      CI_inv = (CI_inv_star + t(CI_inv_star))/2
      logdetCI = logdetCI_star
      PCPt = PCPt_star
      
      l = rep(0,p)
      gamma_l = rep(0,p)
      X_sweep = sweep(X, 2, l, `*`)
      
    }
  }
  
  # sampling step for rho when is != 0 
  if (rho != 0){
    rho_jump = .1
    rho_star = rgamma(1, shape = rho^2/rho_jump^2,
                      rate = rho/rho_jump^2)
    CI_star = (rho_star/rho)*PCPt + diag(n)*sigmay
    chol = chol(CI_star)
    logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
    CI_inv_star = chol2inv(chol)
    
    # TODO: CHECK THIS PLZZ
    loglik = logdetCI_star - logdetCI +
      sum(sapply(1:Ty, function(i) 
        t(mu[ty==i])%*%(CI_inv_star - CI_inv)%*%mu[ty==i]))
    loglik = -0.5*loglik
    
    logprop = dgamma(rho,shape = rho_star^2/rho_jump^2,
                     rate = rho_star/rho_jump^2,log = T) -
      dgamma(rho_star,shape = rho^2/rho_jump^2,
             rate = rho/rho_jump^2,log = T)
    logprior = dgamma(rho_star,1,log = T) -
      dgamma(rho,1,log = T)
    
    logr = loglik + logprop + logprior
    logu = log(runif(1))
    
    if (logr > logu){
      rho = rho_star
      CI = (t(CI_star) + CI_star)/2
      CI_inv = (CI_inv_star + t(CI_inv_star))/2
      logdetCI = logdetCI_star
      PCPt = P%*%CI%*%Pt
    }
  }
  
  params[["rho"]] = rho
  params[["CI"]] = CI
  params[["CI_inv"]] = CI_inv
  params[["logdetCI"]] = logdetCI
  params[["PCPt"]] = PCPt
  
  return(params)
}