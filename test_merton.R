source("MertonProcess.R")
source("plot_simulations.R")

library(tictoc)


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



# par(mfrow=c(3,1))
# 
# plot(t, res$S,pch = 20)
# lines(t, S0*exp((mu-sigma*sigma*0.5)*t), type = 'l', col = 'blue')
# title("S")
# 
# plot(t, res$X, pch = 20)
# lines(t, (mu-sigma*sigma*0.5)*t, type = 'l', col = 'blue')
# title("returns")
# 
# plot(t[-1], cumsum(res$jumps),pch = 20)
# title("Compound poisson")


#######################################################################
# Test speed of mine vs old

plot(t, res$X,pch = 20)


res2 = simulateJump(mu,sigma,lambda,mu_j,sigma_j,TT, t[2]-t[1])
tt = (as.matrix(res2))[,1]
ss = (as.matrix(res2))[,2]
plot(tt,ss,pch = 20)

########################################################################

# simulations
tic("my code")
n = 2000
res =MertonProcess(S0,mu,sigma,lambda,mu_j,sigma_j,TT, t)
t = res$t
sims = res$X

for (i in 2:n) {
  sims = rbind(sims, (MertonProcess(S0,mu,sigma,lambda,mu_j,sigma_j,TT, t))$X)
  print(i)
}

l= length(t)

x11()
hist(sims[,l],breaks = 30)
toc()


tic("other code")
n = 2000
res2 =simulateJump(mu,sigma,lambda,mu_j,sigma_j,TT, t[2]-t[1])
t =(as.matrix(res2))[,1]
sims2 = (as.matrix(res2))[,2]

for (i in 2:n) {
  sims2 = rbind(sims2, (as.matrix(simulateJump(mu,sigma,lambda,mu_j,sigma_j,TT, t[2]-t[1])))[,2])
  print(i)
}

l= length(t)

# plot_simulations(t,sims2)
x11()
hist(sims2[,l],breaks = 30)

toc()

mean(sims[,l])
sd(sims[,l])

mean(sims2[,l])
sd(sims2[,l])


