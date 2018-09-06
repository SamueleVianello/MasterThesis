## Simple Brownian Motion Simulation
mu = 0.5
sigma = 0.5
S0 = 1

dt = 0.1
T = 5

# simulation function
BrM = function(mu,sigma, S0, T, dt){
N = (T/dt)
t = (0:N)*dt
S = rep(0,N+1)
S[1] = S0
for (i in 1:N){
  S[i+1]=S[i]+mu*dt + sigma*sqrt(dt)*rnorm(1)
}
return(rbind(t,S))
}


# simulation
n = 20
res =BrM(mu,sigma,S0,T,dt)
t = res[1,]
sims = res[2,]

for (i in 2:n) {
  sims = rbind(sims, BrM(mu,sigma,S0,T,dt)[2,])
}

# plot 
plot_simulations(t,sims,FALSE)

# check on final distribution
sqrt(var(sims[,dim(sims)[2]]))
mean(sims[,dim(sims)[2]])



##================== Geometric BrM ======================
GBrM = function(mu,sigma, S0, T, dt){
  N = (T/dt)
  t = (0:N)*dt
  S = rep(0,N+1)
  X = rep(0,N+1)
  S[1] = S0
  z = rnorm(N)
  for (i in 2:N){
    X[i+1]=(mu-0.5*sigma^2)*dt + sigma*sqrt(dt)*z[i]
  }
  S=S0*exp(cumsum(X))
  return(rbind(t,S))
}

# simulation
n = 20
res = BrM(mu,sigma,S0,T,dt=0.01)
t = res[1,]
sims = res[2,]

for (i in 2:n) {
  sims = rbind(sims, BrM(mu,sigma,S0,T,dt=0.01)[2,])
}

# plot 
plot_simulations(t,sims,FALSE)


#plot(t,S0*exp((mu-0.5*sigma^2)*t + sigma*sqrt(t)*rnorm(501)),type = 'l')
lines(t,S0*exp((mu-0.5*sigma^2)*t),type = 'l')
