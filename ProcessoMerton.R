ProcessoMerton <-function(drift,sigma,N,lambdaj,mu_xi,h){



  t=(0:N/h)/N
  X=rep(0, N+1)
  F=rep(0, N+1)
  I=rep(0,N)
  X[1]=1
  for(j in 1:N) {
    I[j]=PPgen(h*lambdaj)
    if (I[j]==0){F[j]=0} else {
      F[j]=mu_xi 
    }
    X[j+1]=X[j] + drift*h*X[j]+sigma*sqrt(h)*rnorm(1)*X[j]-F[j]*X[j]
    
  }

return (X)
}