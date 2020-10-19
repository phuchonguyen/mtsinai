source(file.path(getwd(), "sim_001.R"))

# Constants
SEED <- 123
set.seed(SEED)

Ks <- c(1, 2, 5)
Txs <- c(1)
Tys <- c(1)
gammas <- c(TRUE, FALSE)
gammaprobs <- c(0, 1, 0.5)

for(K in Ks) {
  for (t in 1:length(Txs)) {
    Tx <- Txs[t]
    Ty <- Tys[t]
    for (gamma in gammas) {
      if (gamma) {
        for (prob in gammaprobs) {
          sim_001(Tx=1, Ty=1, K=1, gamma=gamma, gamma.prob=prob,
                  niter=20, nburn=10, nthin=5)
        }
      } else {
        sim_001(Tx=1, Ty=1, K=1, gamma=gamma,
                niter=20, nburn=10, nthin=5)
      }
    } 
  }
}









