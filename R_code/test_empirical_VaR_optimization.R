# empirical VaR portfolio optimization

source("MarkowitzMeanVariancePortfolio.R")
source("PortfolioVaROptimization.R")
load("returns.Rda")
load("data.Rda")
load("results.Rda")


N_samples=dim(my_returns)[1]
asset_names = colnames(my_returns[,2*(1:14)])


# rendimenti percentuali:
percentage_returns = exp(my_returns[,2*(1:14)])

hist(percentage_returns$sp500, freq = FALSE, main="Percentage daily S&P500 return")


# ANNUAL percentage simulated return matrix

N_assets = length(w)
N_sim= 1000 # takes almost 45 mins for 10000, 5 min for 1000

sim_mat = matrix(rep(0,N_sim*N_assets), nrow = N_sim)
x=rep(1,N_assets)

t_beg = Sys.time()
for(k in 1:N_sim){
  x = x=rep(1,N_assets)
  samples = sample(N_samples,255, replace = TRUE)
  for( i in 1:255){
    x = x* percentage_returns[samples[i],]
  }
  sim_mat[k,]=as.matrix(x)
}
t_end = Sys.time()

t_end-t_beg

colnames(sim_mat)=asset_names

w= rep(1/N_assets, N_assets)
annual_returns = sim_mat %*% w

hist(log(annual_returns), breaks = 80, main = "Naif allocation annual log-return")
abline(v = quantile(log(annual_returns),probs = 0.05), col = 'blue')

VaR = quantile(1-annual_returns,probs = 1-c(0.01,0.05,0.1))
print(paste("percentage VaR for naif allocation at",  c(0.01,0.05,0.1)*100 , "% level is",VaR))



#################################################################
################### SIMULATING FRONTIER #########################
#################################################################

# Pretty time consuming, a single optimization using 10 repetitions take 2 mins (10000 scenarios)
# So use low number of returns for frontier

returns= c(1.05,1.1,1.15,1.20,1.25)

t_beg = Sys.time()
allocations = matrix(rep(0,N_assets*length(returns)), ncol=N_assets)
VaRs = rep(0,length(returns))
resulting_returns = rep(0,length(returns)) # just as a check
for (i in 1:length(returns)){
  print(returns[i])
  sol_func = OptimalAllocationVaR(simulated_returns = sim_mat, alpha = 0.05, target_return = returns[i], N_rep = 10)
  allocations[i,]=sol_func$allocation
  VaRs[i]=sol_func$objective
  resulting_returns[i]= sol_func$expected_return
}
t_end = Sys.time()
t_end-t_beg


plot(VaRs, returns, type='l')



# single test
sol_func = OptimalAllocationVaR(simulated_returns = sim_mat, alpha = 0.05, target_return = 1.05, N_rep = 10)
sol_func


idx=5
hist(sim_mat%*%allocations[idx,],breaks=100)
abline(v = quantile(sim_mat%*%allocations[idx,],probs = 0.05), col = 'blue')
abline(v = resulting_returns[idx], col='red')
abline(v = mean(sim_mat%*%allocations[idx,]))


ptf_cvar(w= allocations[5,], sims = sim_mat,alpha = 0.05,target_return = 121)


########################## DAILY VaR RETURN ###################





