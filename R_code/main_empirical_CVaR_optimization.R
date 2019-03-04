# empirical VaR portfolio optimization

source("MarkowitzMeanVariancePortfolio.R")
source("PortfolioVaROptimization.R")
load("returns.Rda")
load("data.Rda")
load("results.Rda")


N_samples=dim(my_returns)[1]
N_assets = dim(my_returns)[2]/2 -1 # -1 to exclude vix
asset_names = colnames(my_returns[,2*(1:N_assets)])


percentage_returns = exp(my_returns[,2*(1:N_assets)])

tot_days = dim(percentage_returns)[1]

last_days_to_use = tot_days

daily_return_matrix = as.matrix(percentage_returns[1:last_days_to_use,])


# daily returns to plot frontier and compute allocations
daily_returns = seq(1.00,1.2, by=0.004)^(1/255)
alpha_percentage = 95

daily_returns=sort(daily_returns,decreasing = FALSE)

# including btc
t_beg = Sys.time()
daily_allocations_cvar = matrix(rep(0,N_assets*length(daily_returns)), ncol=N_assets)
colnames(daily_allocations_cvar)=asset_names
daily_CVaRs = rep(0,length(daily_returns))
daily_resulting_returns_cvar = rep(0,length(daily_returns)) # just as a check
for (i in 1:length(daily_returns)){
  print(daily_returns[i]^255)
  sol_func = OptimalAllocationDailyCVaR(daily_return = daily_return_matrix, alpha = 1-alpha_percentage/100, target_return = daily_returns[i], N_rep = 1)
  daily_allocations_cvar[i,]=sol_func$allocation
  daily_CVaRs[i]=sol_func$objective
  daily_resulting_returns_cvar[i]= sol_func$expected_return
}
t_end = Sys.time()
t_end-t_beg



plot(daily_CVaRs, (daily_resulting_returns_cvar)-1, type='l', col = "green", xlab = paste0("Daily CVaR ", alpha_percentage,"%"), ylab = "Daily Returns in %")
title("Efficient CVaR Frontier (daily) including BTC")
grid()


# excluding btc
max_return_no_btc = max(colMeans(percentage_returns[,2:N_assets]))
daily_returns_reduced = daily_returns[which(daily_returns <= max_return_no_btc)]

t_beg = Sys.time()
daily_allocations_cvar_no_btc = matrix(rep(0,N_assets*length(daily_returns_reduced)), ncol=N_assets)
colnames(daily_allocations_cvar_no_btc) = asset_names
daily_CVaRs_no_btc = rep(0,length(daily_returns_reduced))
daily_resulting_returns_cvar_no_btc = rep(0,length(daily_returns_reduced)) # just as a check
for (i in 1:length(daily_returns_reduced)){
  print(daily_returns_reduced[i]^255)
  sol_func = OptimalAllocationDailyCVaR(daily_return = daily_return_matrix[,2:N_assets], alpha = 1-alpha_percentage/100, target_return = daily_returns_reduced[i], N_rep = 1)
  daily_allocations_cvar_no_btc[i,2:N_assets]=sol_func$allocation
  daily_CVaRs_no_btc[i]=sol_func$objective
  daily_resulting_returns_cvar_no_btc[i]= sol_func$expected_return
}
t_end = Sys.time()
t_end-t_beg

plot(daily_CVaRs_no_btc, (daily_resulting_returns_cvar_no_btc)-1, type='l', col = "orange", xlab = paste0("Daily CVaR ", alpha_percentage,"%"), ylab = "Daily Returns in %")
title("Efficient CVaR Frontier (daily) excluding BTC")
grid()



# polish allocation data
daily_allocations_cvar[which(abs(daily_allocations_cvar)<1e-8)] =0
daily_allocations_cvar_no_btc[which(abs(daily_allocations_cvar_no_btc)<1e-10)] =0 

# aggregate results to be saved as csv file
res_btc = cbind(daily_resulting_returns_cvar-1,daily_CVaRs,daily_allocations_cvar)
colnames(res_btc)[c(1,2)] = c("return_daily", "cvar_daily")
res_no_btc = cbind(daily_resulting_returns_cvar_no_btc-1,daily_CVaRs_no_btc,daily_allocations_cvar_no_btc )
colnames(res_no_btc)[c(1,2)] = c("return_daily", "cvar_daily")


# # save to file
# write.csv(file = paste0("allocation_cvar",alpha_percentage, ".csv"), x = res_btc)
# write.csv(file = paste0("allocation_cvar",alpha_percentage, "_no_btc.csv"), x = res_no_btc)
