MertonProcess <-function(mu,sigma,lambda,mu_j,sigma_j,TT, dt){
  # INPUT
  # mu:
  # sigma:
  # lambda: Poisson parameter
  # mu_j: mean of (log of) jumps
  # sigma_j: volatility of (log of) jumps

  # VARIABLES
  # X: log of 

  t=(0:N/h)/N
  X=rep(0, N+1)
  F=rep(0, N+1)
  I=rep(0,N)
  X[1]=1
  for(j in 1:N) {
    I[j]=PPgen(h*lambda)
    if (I[j]==0){F[j]=0} else {
      F[j]=mu_j 
    }
    X[j+1]=X[j] + mu*h*X[j]+sigma*sqrt(h)*rnorm(1)*X[j]-F[j]*X[j]
    
  }

return (X)
}