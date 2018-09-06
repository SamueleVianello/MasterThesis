CompoundPoissonProcess=function(lambda, TT, mu_j, sigma_j, t){
  # Simulates a compound poisson process using alg 6.2
  # in Cont Tankov (2004). Jump size is assumed normally
  # distributed N(mu_j,sigma_j^2). Assumed NO drift in between jumps.
  
  # INPUT:
  # lambda: poisson parameter
  # TT: final time in simulation interval [0,TT]
  # mu_j: mean of jump size
  # sigma_j: vol of jump size
  # t: time instants for simulation
  
  # OUTPUT
  # X: compound poisson process valued at t
  
  # total number of jumps in [0,TT]
  N = rpois(1,lambda*TT) 
  # jump times
  t_j = sort(runif(N)*TT)
  # jump sizes
  J = rnorm(N,mean = mu_j, sd = sigma_j)
  
  # simulation
  jump_sum = 0;
  N_j = 1;
  l = length(t)
  X = rep(0,l)
  
  for (i in 1:l){
    if (N_j <=N & t[i] > t_j[N_j]){
      jump_sum = jump_sum + J[N_j]
      N_j = N_j + 1
    }
    X[i] = jump_sum
  }

  return(X)
}