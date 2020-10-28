############################################
# Simulation from paper                    #
############################################
n <- 50  # number of samples
p <- 10  # number of predictors
q <- 3   # number of outcomes
# Active number of predictors is 2
p.act <- 2

#############################
# Generate design matrix X. #
#############################
times <- 1:p
rho <- 0.85
H <- abs(outer(times, times, "-"))
V <- rho^H
mu <- rep(0, p)
# Rows of X are simulated from MVN(0,V)
X <- MASS::mvrnorm(n, mu, V)
# Center X
X <- scale(X, center=TRUE, scale=FALSE)


########################################
# Generate true error covariance Sigma #
########################################
sigma.sq=2
times <- 1:q
H <- abs(outer(times, times, "-"))
Sigma <- sigma.sq * rho^H

###########################
# Generate noise matrix E #
###########################
mu <- rep(0,q)
E <- MASS::mvrnorm(n, mu, Sigma)

#########################################
# Generate true coefficient matrix B_0. #
#########################################
# Entries in nonzero rows are drawn from Unif[(-5,-0.5)U(0.5,5)]
# Option 1: B[i,j] iid from U(-5, 4)
#B.act <- runif(p.act*q,-5,4)
# Option 2: B[i,] iid from N(m, Sigma), E[i,] ~ N(0, Sigma)
B.act <- MASS::mvrnorm(p.act, seq(-5, 4, length.out = q), Sigma)  # Correlated B
disjoint <- function(x){
  if(x <= -0.5) return(x)
  else return(x+1)
}
B.act <- matrix(sapply(B.act, disjoint),p.act,q)
# Set rest of the rows equal to 0
B.true <- rbind(B.act,matrix(0,p-p.act,q))
B.true <- B.true[sample(1:p),] # permute the rows

# Option 3: B from N(M, V, Sigma), Var(X) = V, Var(E) = Sigma
# B.act <- MBSP::matrix.normal(M = matrix(sapply(runif(p*q,-5,4), disjoint), p, q),
#                              U = V, V = Sigma)
# B.act[(p.act+1):p, ] <- 0
# B.true <- B.act[sample(1:p),]


##############################
# Generate response matrix Y #
##############################
Y <- crossprod(t(X),B.true) + E


#########################################
# Run the MBSP model on synthetic data. #
#########################################
# For optimal estimation, change max.steps to 15,000
# and change burnin to be 5000
mbsp.model = MBSP::mbsp.tpbn(X=X, Y=Y, max.steps = 15000, burnin=5000)

par(mfrow=c(1,2))
image(mbsp.model$B.est)
image(B.true)


###############################################################
# Muri model, essentially the same as MBSP but with HIW prior #
#             also option to select each variable available   #
###############################################################
data<- list(X=X, Y=Y)
# Shrinkage prior
samps <- Muri(15000, data)
# Spike and slab prior
samps.select <- Muri(15000, data, method="select")


#################################################
#  Some summary                                 #
#################################################
thin <- seq(5000, 15000, 10)
plot(samps$B[thin,2,1], type="l")
plot(samps$Sigma[thin,1,1], type="l")
B.post <- apply(samps$B[thin,,], c(2,3), mean)
par(mfrow=c(1,2))
image(B.post)
image(B.true)

plot(samps.select$B[thin,2,1], type="l")
plot(samps.select$Sigma[thin,1,1], type="l")
B.post <- apply(samps.select$B[thin,,], c(2,3), median)
par(mfrow=c(1,2))
image(B.post)
image(B.true)

Sigma.post <- apply(samps$Sigma[thin,,], c(2,3), mean)
par(mfrow=c(1,3))
image(cov(Y))
image(Sigma.post)
image(Sigma)

par(mfrow=c(1,3))
image(cov(B.true))
image(cov(mbsp.model$B.est))
image(cov(B.post))