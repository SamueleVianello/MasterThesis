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

returns= c(1.05,1.1,1.15,1.20,1.25,1.30,1.35,1.40)



#--------------------------------------------------------------
##################### ANNUAL VaR/CVaR RETURN ##################
#--------------------------------------------------------------
## VAR

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


plot(VaRs, resulting_returns, type='l')

# lines(VaRs, log(resulting_returns), type='l', col='blue')

# single test
sol_func = OptimalAllocationVaR(simulated_returns = sim_mat, alpha = 0.05, target_return = 1.05, N_rep = 10)
sol_func


idx=5
hist(sim_mat%*%allocations[idx,],breaks=100)
abline(v = quantile(sim_mat%*%allocations[idx,],probs = 0.05), col = 'blue')
abline(v = resulting_returns[idx], col='red')
abline(v = mean(sim_mat%*%allocations[idx,]))




# CVaR

t_beg = Sys.time()
allocations_cvar = matrix(rep(0,N_assets*length(returns)), ncol=N_assets)
CVaRs = rep(0,length(returns))
resulting_returns_cvar = rep(0,length(returns)) # just as a check
for (i in 1:length(returns)){
  print(returns[i])
  sol_func = OptimalAllocationVaR(simulated_returns = sim_mat, alpha = 0.05, target_return = returns[i], N_rep = 10, CVaR = TRUE)
  allocations_cvar[i,]=sol_func$allocation
  CVaRs[i]=sol_func$objective
  resulting_returns_cvar[i]= sol_func$expected_return
}
t_end = Sys.time()
t_end-t_beg

plot(CVaRs, log(resulting_returns_cvar), type='l')






#--------------------------------------------------------------
####################### DAILY VaR/CVaR RETURN #################
#--------------------------------------------------------------
tot_days = dim(percentage_returns)[1]
last_5y = 255*5
last_5y_returns = as.matrix(percentage_returns[1:last_5y,])

sol_func = OptimalAllocationDailyVaR(daily_return =last_5y_returns, alpha = 0.05, target_return = 1.2^(1/255), N_rep = 10)
sol_func


w = rep(1/N_assets,N_assets)
# w = c(0, rep(1/(N_assets-1),N_assets-1))
daily_ptf_cvar(w = w, daily_return = last_5y_returns, alpha = 0.05,target_return = 11)
sum(w*colMeans(last_5y_returns))


# daily returns to plot frontier
daily_returns = returns^(1/255)


t_beg = Sys.time()
daily_allocations_cvar = matrix(rep(0,N_assets*length(returns)), ncol=N_assets)
daily_CVaRs = rep(0,length(returns))
daily_resulting_returns_cvar = rep(0,length(returns)) # just as a check
for (i in 1:length(returns)){
  print(returns[i])
  sol_func = OptimalAllocationDailyVaR(daily_return = last_5y_returns, alpha = 0.05, target_return = daily_returns[i], N_rep = 10)
  daily_allocations_cvar[i,]=sol_func$allocation
  daily_CVaRs[i]=sol_func$objective
  daily_resulting_returns_cvar[i]= sol_func$expected_return
}
t_end = Sys.time()
t_end-t_beg

colnames(daily_allocations_cvar)=asset_names

plot(daily_CVaRs, (daily_resulting_returns_cvar), type='l')


#------------------------------------------------------------------
################### PLOTTING PORTFOLIO COMPOSITION ################
#------------------------------------------------------------------
# ggplot2 library
library(ggplot2)
library(RColorBrewer)


l=length(returns)


# annual CVaR from empirical simulation
Asset = rep(asset_names, l)
Return = rep(resulting_returns_cvar, each = 14)
Values = drop(matrix(t(allocations_cvar), nrow = 1, byrow = TRUE))

# colorRampPalette(brewer.pal(9, "Spectral"))(14)

data_annual <- data.frame(Asset,Return,Values)
ggplot(data_annual, aes(x=Return, y=Values, fill=Asset)) + 
  geom_area(alpha=1 , size=1, colour="black") +
  ggtitle("Empirical CVaR allocation")+
  scale_fill_manual(values =colorRampPalette(brewer.pal(9, "Paired"))(14) )+
  scale_y_continuous(breaks=seq(0,1,by=0.1))




# Daily CVaR from previous 5 years

Asset_daily = rep(asset_names, l)
Return_daily = rep(daily_resulting_returns_cvar, each = 14)
Values_daily = drop(matrix(t(daily_allocations_cvar), nrow = 1, byrow = TRUE))



# colorRampPalette(brewer.pal(9, "Spectral"))(14)

data_daily <- data.frame(Asset_daily,Return_daily,Values_daily)
ggplot(data_daily, aes(x=Return, y=Values, fill=Asset)) + 
  geom_area(alpha=1 , size=1, colour="black") +
  ggtitle("Daily CVaR allocation")+
  scale_fill_manual(values =colorRampPalette(brewer.pal(9, "Paired"))(14) )+
  scale_y_continuous(breaks=seq(0,1,by=0.1))+
  scale_x_continuous(sec.axis = ~.^(255))

