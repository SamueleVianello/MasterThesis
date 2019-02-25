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

# ***** LOG-RETURNS ********
# SS = cov(my_returns[1:N_samples,2*(1:N_assets)]) * 255
# expected_return_sample = colMeans(my_returns[1:N_samples,2*(1:N_assets)]) * 255
# y_lim = c(-0.05,0.4)
# max_r = 0.4

# ***** PERCENTAGE RETURNS ******
expected_return_sample = colMeans(exp(my_returns[1:N_samples,2*(1:N_assets)]))^255
SS = cov(exp(my_returns[1:N_samples,2*(1:N_assets)]))*255
max_r = 1.6
y_lim = c(1.00,max_r)




# Plots the efficient frontiers, w/ and w/o Bitcoin, w/ and w/o shortselling
eff_front = PlotEfficientFrontier(expected_return_sample, SS, min_r = 1, max_r=max_r,
                                  exclude_btc = TRUE, add_no_short_sell = T, full_plot = FALSE)


# dev.copy2pdf(file = "efficient_frontier.pdf", height = 7, width=7 )
# dev.off()



# # target return
# target = 0.2
# 
# w = OptimalAllocation(r=expected_return_sample,S=SS, target_return = target)
# w_no_btc = OptimalAllocation(r=expected_return_sample[2:14],S=SS[2:14,2:14], target_return = target)
# 
# cbind(w, c(0,w_no_btc))

# percentage ****************
targets = seq(1.025, 1.17, by = 0.005)

# log-returns *******************
# targets = seq(0.05,0.30, by=0.01)


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

# plot(targets,alloc_btc[1,], type='l', col=rainbow(14)[1], ylim = c(-1,1))
# for(i in 2:14){
#   lines(targets, alloc_btc[i,], col= rainbow(14)[i])
# }


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
