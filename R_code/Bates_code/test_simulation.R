source("SimulationBates.R")
load("returns.Rda")

params = c(0.127 , 15.67, 0.01865, 0.764, -0.00042)
check_feller(params = params)


dt = 1/255
TT = 1

Nsim = 2000

set.seed(1234)
sims = simulate_mv_heston(Nasset = 1,mu = 0.127,k = 15.67, eta = 0.01865,theta = 0.764,rho = -0.00042,
                          S0 = 1, V0 = 0.017, dt = dt,final_t = TT,Nsim = Nsim)

sims$return
sims$variance

tt = (0:(TT/dt))*dt

plot(tt, sims$returns, type = 'l')

hist(sims$returns[,256])
mean(sims$returns[,256]) + var(sims$returns[,256])*0.5
var(sims$returns[,256])

matplot(tt, t(sims$variance$V_1[1:2,]), type = 'l')



mvsim=simulate_mv_heston(3,0.127*c(1,0.5,0.2), 15.67*c(1,0.5,0.2), 0.01865*c(1,0.5,0.2), 0.764*c(1,0.5,0.2), -0.003*c(1,0.5,0.2), 
                         matrix(c(1,0.3,-0.4,0.3,1,0.7,-0.4,0.7,1),ncol = 3,nrow = 3), 1*c(1,1,1), 0.02*c(1,0.5,0.2), dt=1/255, final_t =1, Nsim = 200 )
correlation_MC_estimation(mvsim)


matplot(tt, cbind(mvsim$return$x_1[1,],mvsim$return$x_2[1,],mvsim$return$x_3[1,]), type = 'l', main = "Assets log-returns dynamics")
matplot(tt, cbind(mvsim$variance$V_1[1,],mvsim$variance$V_2[1,],mvsim$variance$V_3[1,]), type = 'l',  main = "Assets variance dynamics")





TT =1
mvsim=simulate_mv_bates(2,0.127*c(1,0.5), 1.67*c(1,0.5), 0.01865*c(1,0.5), 0.1764*c(1,0.5), -0.003*c(1,0.5),
                        -0.1*c(1,0.5), c(0.05,0.05), c(0, 0),
                         matrix(c(1,-0.3584472,-0.3584472,1),ncol = 2,nrow = 2), 
                        1*c(1,1), 0.02*c(1,0.5), 
                        dt =1/255, final_t = TT, 
                        Nsim = 10000 )
correlation_MC_estimation(mvsim)

mean(mvsim$return$x_1[,255*TT])
mean(mvsim$return$x_1[,255*TT]) + var(mvsim$return$x_1[,255*TT])*0.5

var(mvsim$return$x_1[,255*TT])

mean(mvsim$variance$V_1[,255*TT])
hist(mvsim$return$x_1[,255*TT])

tt = (0:(TT/dt))*dt
matplot(tt, exp(t(mvsim$return$x_1[1:100,])), type = 'l', main = "Assets log-returns dynamics")
matplot(tt, cbind(mvsim$variance$V_1[1,],mvsim$variance$V_2[1,]), type = 'l',  main = "Assets variance dynamics")
abline(h=0)



######### 
xx = seq(-1,1,length.out = 200)
yy = -xx + log(1+xx)
plot(xx,yy, type='l', ylim=c(-1,0))
grid()




# testing expected model correlation estimation
TT=1

model_rho = calibrate_correlation(sample_corr = -0.03, 2,mu=0.127*c(1,0.5), k=1.67*c(1,0.5), eta=0.01865*c(1,0.5), theta=0.1764*c(1,0.5),rho= -0.003*c(1,0.5),
  mu_j=-0.1*c(1,0.5), sigma_j= c(0.05,0.05), lambda=c(0.5, 0.2),
  V0= 0.02*c(1,0.5), 
  dt =1/255, final_t = TT, 
  Nsim = 2000 )

model_rho


mvsim=simulate_mv_bates(2,0.127*c(1,0.5), 1.67*c(1,0.5), 0.01865*c(1,0.5), 0.1764*c(1,0.5), -0.003*c(1,0.5),
                        -0.1*c(1,0.5), c(0.05,0.05), c(0.5, 0.2),
                        matrix(c(1,model_rho,model_rho,1),ncol = 2,nrow = 2), 
                        1*c(1,1), 0.02*c(1,0.5), 
                        dt =1/255, final_t = TT, 
                        Nsim = 10000 )
correlation_MC_estimation(mvsim)



# test full matrix correlation estimation -----------------------------------

# asset parameters
par1 = c(0.000673,19.7,0.0294,0.8196,-0.422)   #bric
par2 = c(0.127,17.07,0.01858, 0.7962,-0.00076) #sp500
par3 = c(0.04756,1.22,0.03807,0.265, -0.5826)  #eurostoxx

pars = rbind(par1,par2,par3)

mu=pars[,1]
k =pars[,2]
eta=pars[,3]
theta=pars[,4]
rho=pars[,5]

# sample corr from given returns (can be substituted to any other choice)
sample_corr = cor(cbind(my_returns[,c(4,6,8)]))

# estimating the correlation matrix for the model
model_c = calibrate_full_correlation(sample_corr_matrix = sample_corr,Nasset = 3,
                                     mu = mu,k = k,eta = eta,theta = theta,rho = rho,
                                     V0 = eta*1.1, dt =1/255, final_t = 1,Nsim = 1000 )
model_c

# Simulating assets to get the time series of returns
sim_3 = mvsim=simulate_mv_heston(Nasset = 3,mu = mu,k = k,eta = eta,theta = theta,rho = rho, CorrMatrix = model_c,
                                S0 =1*c(1,1,1), V0 = eta*1.1, dt=1/255, final_t =1, Nsim = 200 )

# Comparing result to sample correlation
correlation_MC_estimation(mvsim)
sample_corr


# 
# 
# ############ testing ESGtoolkit ################
# 
# library(ESGtoolkit)
# 
# # Spot variance
# V0 <- 0.022
# # mean-reversion speed
# kappa <- 15.67
# # long-term variance
# theta <- 0.01865
# # volatility of volatility
# volvol <- 0.764
# # Correlation between stoch. vol and prices
# rho <- -0.00042
# 
#  
# # # Intensity of the Poisson process
# # lambda <- 0.3635
# # # mean and vol of the merton jumps diffusion
# # mu.J <- -0.2459
# # sigma.J <- 0.2547/100
# # m <- exp(mu.J + 0.5 * (sigma.J^2)) - 1
# 
# # Initial stock price
# S0 <- 100
# # Initial short rate
# r0 <- 0.127
# 
# n <- 1000
# horizon <- 1
# 
# freq <- "weekly"
# # Simulation of shocks, with antithetic variates
# shocks <- simshocks(n = n, horizon = horizon, frequency = freq, method = "anti",
#                     family = 1, par = rho)
# # Vol simulation
# sim.vol <- simdiff(n = n, horizon = horizon, frequency = freq, model = "CIR",
#                    x0 = V0, theta1 = kappa * theta, theta2 = kappa, theta3 = volvol, eps = shocks[[1]])
# # Plotting the volatility (only for a low number of simulations)
# esgplotts(sim.vol)
# 
# 
# # prices simulation
# sim.price <- simdiff(n = n, horizon = horizon, frequency = freq, model = "GBM",
#                      x0 = 1, theta1 = r0 - lambda * m, theta2 = sim.vol, 
#                      # lambda = lambda, mu.z = mu.J, sigma.z = sigma.J, 
#                      eps = shocks[[2]])
# 
# final_ret = log(sim.price[53,])
# hist(final_ret)
# mean(final_ret) + 0.5*var(final_ret)
# var(final_ret)
# 
# 
# x11()
# par(mfrow = c(2, 1))
# matplot(time(sim.price), sim.price, type = "l", main = "with matplot")
# esgplotbands(sim.price, main = "with esgplotbands", xlab = "time", ylab = "values")
