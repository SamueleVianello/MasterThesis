source("MertonProcess.R")
source("plot_simulations.R")


lambda = 10
TT = 1
mu=1
sigma = 0.31
mu_j =0
sigma_j = 0.1
Nsim = 1000
t = (0:(TT*Nsim))/(Nsim)
S0 = 100


# plot(t,CompoundPoissonProcess(lambda,TT, mu_j, sigma_j,t))



res =MertonProcess(S0,mu,sigma,lambda,mu_j,sigma_j,TT, t)



par(mfrow=c(3,1))

plot(t, res$S,pch = 20)
lines(t, S0*exp((mu-sigma*sigma*0.5)*t), type = 'l', col = 'blue')
title("S")

plot(t, res$X, pch = 20)
lines(t, (mu-sigma*sigma*0.5)*t, type = 'l', col = 'blue')
title("returns")

plot(t[-1], cumsum(res$jumps),pch = 20)
title("Compound poisson")


