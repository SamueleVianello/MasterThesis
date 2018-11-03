source("MarkowitzMeanVariancePortfolio.R")
load("returns.Rda")
load("data.Rda")
load("results.Rda")



##########################################
##### test on ourcalibrated assets #######
##########################################
library(pracma)


# attach(my_returns)
source("MarkowitzMeanVariancePortfolio.R")

# Number of samples to consider from latest {max is dim(my_returns)[1]}
# N_samples = 2*255
N_samples = dim(my_returns)[1]
N_assets = dim(my_returns)[2]/2

# From sample
SS = cov(my_returns[1:N_samples,2*(1:N_assets)]) * 255
# expected_return_sample = colMeans(my_returns[1:N_samples,2*(1:N_assets)]) * 255
# y_lim = c(-0.05,0.4)
# max_r = 0.4

# returns as percentage:
# expected_return_sample = colMeans(exp(my_returns[1:N_samples,2*(1:N_assets)]))^255
# y_lim = c(1.00,1.4)
# max_r = 1.4

# yearly from model
# expected_return =  results$parameters[1,]
# y_lim = c(-0.05,0.4)
# max_r = 0.4
SS = results$covariance

# percentage return from model
expected_return = exp(results$parameters[1,])
y_lim = c(1.00,1.4)
max_r = 1.4

# unconstrained = short sales are allowed for every asset
res1 = EfficientFrontier(expected_return_sample,SS, max_r = max_r)
res2 = EfficientFrontier(expected_return_sample[2:N_assets],SS[2:N_assets,2:N_assets], max_r = max_r)


x11()

#windows(width = 10,height = 8)
plot(res1$sigma,res1$expected_return, type = 'l',col='darkgreen', ylim = y_lim, xlab = "Volatility", ylab = "Returns")
lines(res2$sigma,res2$expected_return, col = 'red')
points(sqrt(diag(SS)),expected_return_sample, pch='+', col = 'blue')
text(sqrt(diag(SS)),expected_return_sample, labels = colnames(my_returns[,2*(1:N_assets)]),pos = 3)
grid()



# constrained = no short sales for given assets
res_btc = EfficientFrontier(r=expected_return_sample, S=SS, full = FALSE, N=100, no_short_sales = 1:N_assets, max_r = max_r)
res_no_btc = EfficientFrontier(r= expected_return_sample[2:N_assets], S=SS[2:N_assets,2:N_assets],full = FALSE, N =200, no_short_sales = 1:(N_assets-1))

lines(res_btc$sigma, res_btc$expected_return, col= 'green')
lines(res_no_btc$sigma[2:201], res_no_btc$expected_return[2:201], col = 'orange')

title(main = "Efficient Markowitz Mean Variance Frontier")
legend("topleft", legend = c("btc", "NO btc", "btc, NO short sale","NO btc, NO short sale"),
       col=c("darkgreen","red","green","orange"), lwd = 3, lty = c(1,1,1,1), cex=0.75)


# dev.copy2pdf(file = "chiappedacciaio_logreturn.pdf")
# dev.off()


# # target return
# target = 0.2
# 
# w = OptimalAllocation(r=expected_return_sample,S=SS, target_return = target)
# w_no_btc = OptimalAllocation(r=expected_return_sample[2:14],S=SS[2:14,2:14], target_return = target)
# 
# cbind(w, c(0,w_no_btc))

# percentage ****************
targets = seq(1.05,1.3, by = 0.025)

# log-returns *******************
# targets = seq(0.05,0.30, by=0.01)


l=length(targets)
alloc_btc = zeros(N_assets,l)
rownames(alloc_btc)=colnames(my_returns[,2*(1:N_assets)])
alloc_no_btc = zeros(N_assets-1,l)
rownames(alloc_no_btc)=colnames(my_returns[,2*(N_assets)])

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

