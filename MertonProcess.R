MertonProcess <-function(S0,mu,sigma,lambda,mu_j,sigma_j,TT, t){
  # Simulates a trajectory of the merton process with 
  # given parameters on a given time grid.
  # Jump size is assumed to be normally distributed.
  
  # INPUT
  # S0: initial value of the asset
  # mu: drift of continuous component 
  # sigma: volatility of continuous component
  # lambda: Poisson parameter
  # mu_j: mean of (log of) jumps
  # sigma_j: volatility of (log of) jumps

  
  # OUTPUT
  # S: Merton process trajectory 
  # X: log-returns of process
  
  source('CompoundPoissonProcess.R')
  
  # time increments
  dt = diff(t) 
  # standard gaussians
  l = length(t)
  dz = rnorm(l-1)
  # jump increments
  dj = diff(CompoundPoissonProcess(lambda = lambda, TT = TT, mu_j=mu_j, sigma_j = sigma_j, t = t))
      # _________________________________________________________________________________
      # it could be a better/faster way to simulate jump increments as poisson(lambda*dt)
  
  
  # simulation of log-returns
  
  X=rep(0,l)
  for (i in 1:(l-1)){
    X[i+1]= X[i] + (mu-sigma*sigma*0.5)*dt[i] +sigma*sqrt(dt[i])*dz[i] + dj[i]
  }
  
  S = S0*exp((X))
  
  return (list(t = t, S = S, X = X, jumps = dj))
}