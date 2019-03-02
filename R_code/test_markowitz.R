source("MarkowitzMeanVariancePortfolio.R")
load("returns.Rda")
load("data.Rda")
load("results.Rda")



##########################################
##### compute efficient frontier #########
##########################################
library(pracma)


# attach(my_returns)
source("MarkowitzMeanVariancePortfolio.R")

# Number of samples to consider from latest {max is dim(my_returns)[1]}
N_samples = dim(my_returns)[1]
# Number of assets
N_assets = dim(my_returns)[2]/2 -1 # -1 to exclude VIX from our analysis


# ***** PERCENTAGE RETURNS ******
expected_return_sample = colMeans(exp(my_returns[1:N_samples,2*(1:N_assets)]))^255
SS = cov(exp(my_returns[1:N_samples,2*(1:N_assets)]))*255
max_r = 1.6
y_lim = c(1.00,max_r)




# Plots the efficient frontiers, w/ and w/o Bitcoin, w/ and w/o shortselling
eff_front = PlotEfficientFrontier(expected_return_sample, SS, min_r = 1, max_r=max_r,
                                  exclude_btc = T, add_no_short_sell = T, full_plot = F)


# dev.copy2pdf(file = "efficient_frontier.pdf", height = 7, width=7 )
# dev.off()



##################### ALLOCATIONS FOR GIVEN VOLATILITY ###################################

vol_target = seq(from = 0.0275,to = 0.12, by = 0.0025)


allocations_vol = zeros(length(vol_target), N_assets) 
colnames(allocations_vol) = colnames(my_returns[,2*(1:N_assets)])
returns_vol = zeros(length(vol_target),1)
allocations_vol_no_btc = zeros(length(vol_target), N_assets) 
colnames(allocations_vol_no_btc) = colnames(my_returns[,2*(1:N_assets)])
returns_vol_no_btc = zeros(length(vol_target),1)

for (i in 1:length(vol_target)) {
  allocations_vol[i,] = OptimalAllocation(r=expected_return_sample,S=SS, sd=vol_target[i], no_short_sales=1:N_assets)
  returns_vol[i,1] = sum(allocations_vol[i,]*expected_return_sample)

  allocations_vol_no_btc[i,2:N_assets] = OptimalAllocation(r=expected_return_sample[2:N_assets],S=SS[2:N_assets,2:N_assets], sd=vol_target[i], no_short_sales=1:(N_assets-1))
  returns_vol_no_btc[i,1] = sum(allocations_vol_no_btc[i,]*expected_return_sample)
}

# polish data for small allocation ( allocation of 1e-10 set to zero)
allocations_vol[which(abs(allocations_vol)<1e-10)] = 0
allocations_vol_no_btc[which(abs(allocations_vol_no_btc)<1e-10)] = 0

points(vol_target,returns_vol_no_btc-1,col="orange")
points(vol_target,returns_vol-1,col="green")

res_btc = cbind(returns_vol-1,vol_target, allocations_vol)
colnames(res_btc)[c(1,2)] = c("exp_return","volatility")
res_no_btc = cbind(returns_vol_no_btc-1,vol_target, allocations_vol_no_btc)
colnames(res_no_btc)[c(1,2)] = c("exp_return","volatility")

write.csv(file = "allocation_on_vol.csv", x = res_btc)
write.csv(file = "allocation_on_vol_no_btc.csv", x = res_no_btc)



##################### ALLOCATIONS FOR GIVEN RETURNS #####################################

# percentage ****************
targets = seq(1.025, 1.17, by = 0.005)



l=length(targets)
alloc_btc = zeros(N_assets,l)
rownames(alloc_btc)=colnames(my_returns[,2*(1:N_assets)])
alloc_no_btc = zeros(N_assets-1,l)
rownames(alloc_no_btc)=colnames(my_returns[,2*(2:N_assets)])

for(i in 1:l){
  # print(paste(i, targets[i]))
  alloc_btc[,i]= OptimalAllocation(r=expected_return_sample,S=SS, target_return = targets[i], no_short_sales = 1:N_assets)
  alloc_no_btc[,i]= OptimalAllocation(r=expected_return_sample[2:N_assets], S=SS[2:N_assets,2:N_assets], target_return = targets[i], no_short_sales = 1:(N_assets-1))
}

alloc_btc
alloc_no_btc

# ggplot2 library
library(ggplot2)
library(RColorBrewer)




Asset = rep(rownames(alloc_btc), l)
Return = rep(targets, each = N_assets)
Values = drop(matrix(alloc_btc, nrow = 1))




colorRampPalette(brewer.pal(9, "Spectral"))(N_assets)

data <- data.frame(Asset,Return,Values)
ggplot(data, aes(x=Return, y=Values, fill=Asset)) + 
  geom_area(alpha=1 , size=1, colour="black") +
  ggtitle("Markowitz Optimal Allocation Without Shortsellling")+
  scale_fill_manual(values =colorRampPalette(brewer.pal(9, "Paired"))(N_assets) )


volatilities_alloc = rep(0,length(targets)) 

for(i in 1:length(targets)){
  volatilities_alloc[i]= sqrt(t(alloc_btc[,i]) %*% SS %*% alloc_btc[,i])
}



##################################################
with_btc = rbind(rep(0,length(alloc_no_btc)), alloc_no_btc)

with_btc[which(with_btc<1e-6)]=0

Asset = rep(rownames(alloc_btc), l)
Return = rep(targets, each = N_assets)
Values = drop(matrix(with_btc, nrow = 1))




colorRampPalette(brewer.pal(9, "Spectral"))(N_assets)

data <- data.frame(Asset,Return,Values)
ggplot(data, aes(x=Return, y=Values, fill=Asset)) + 
  geom_area(alpha=1 , size=1, colour="black") +
  ggtitle("Markowitz Optimal Allocation Without Bitcoin")+
  scale_fill_manual(values =colorRampPalette(brewer.pal(9, "Paired"))(N_assets) )


volatilities_alloc_no_btc = rep(0,length(targets)) 

for(i in 1:length(targets)){
  volatilities_alloc_no_btc[i]= sqrt(t(alloc_no_btc[,i]) %*% SS[2:N_assets,2:N_assets] %*% alloc_no_btc[,i])
}
