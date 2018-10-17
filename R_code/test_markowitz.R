source("MarkowitzMeanVariancePortfolio.R")
load("returns.Rda")
load("data.Rda")
load("results.Rda")


r =matrix( c(0.03,0.05,0.08),ncol = 1)# vector of returns
sd= diag(c(0.12,0.1,0.2))
corr= matrix(c(   1,-0.1, 0.4,
               -0.1,   1, 0.3,
                0.4, 0.3,   1),3,3)
S = sd%*%corr%*%sd# covariance matrix
invS = solve(S)
e = matrix(rep(1,length(r)),nrow = length(r),ncol = 1) # unit vector

expected_return = 0.4

a = drop(t(e)%*% invS %*% e)
b = drop(t(e)%*% invS %*% r)
c = drop(t(r)%*% invS %*% r)
d = drop(a*c - b^2)

mu = drop((a*expected_return - b)/d)
lambda = drop( (d - a*b*expected_return + b^2)/(a*d))

w_opt = invS %*% (e*lambda + mu*r)

sigma_opt = sqrt(t(w_opt)%*%S%*%w_opt)
sigma = sqrt( (a*expected_return^2 - 2*b*expected_return + c)/d)
sigma
sigma_opt

xx= seq(from = 0,to = 0.1,length.out = 100)
yy = sqrt( (a*xx^2 - 2*b*xx + c)/d)

plot(yy,xx,type ='l')
points(diag(sd),r,col='blue')


EfficientFrontier(r,S,full = TRUE, plot = TRUE)


# test with no short sales

library(quadprog)

D = 2*S
d = matrix(rep(0,length(r)),ncol = 1)

A = rbind(t(r),rep(1,length(r)),
          diag(length(r)))
A
b= c(expected_return, 1, rep(0,length(r)))

sol = solve.QP(Dmat = D, dvec = (d), Amat = t(A), bvec = t(b), meq = 2)



res = EfficientFrontier_constr(r,S, plot=FALSE)
res$sigma
res$expected_return

lines(res$sigma,res$expected_return, type = 'l', col='red')


##########################################
##### test on ourcalibrated assets #######
##########################################
library(pracma)

attach(my_returns)
source("MarkowitzMeanVariancePortfolio.R")

# From sample
expected_return_sample = colMeans(my_returns[,2*(1:14)]) * 255
SS = cov(my_returns[,2*(1:14)]) * 255

# yearly from model
# expected_return =  colSums(results$full_mu)/6 + colSums(results$full_theta)/6 * colSums(results$full_lambda)/6
# SS = results$covariance

# unconstrained = short sales are allowed for every asset
res1 = EfficientFrontier(expected_return_sample,SS, max_r = 0.4)
res2 = EfficientFrontier(expected_return_sample[2:14],SS[2:14,2:14], max_r = 0.4)


#x11()


windows(width = 10,height = 8)
plot(res1$sigma,res1$expected_return, type = 'l',col='darkgreen', ylim = c(-0.1, 0.4), xlab = "Volatility", ylab = "Log-Returns")
lines(res2$sigma,res2$expected_return, col = 'red')
points(sqrt(diag(SS)),expected_return_sample, pch='+', col = 'blue')
text(sqrt(diag(SS)),expected_return_sample, labels = colnames(my_returns[,2*(1:14)]),pos = 3)
grid()



# constrained = no short sales for given assets
res_btc = EfficientFrontier(r=expected_return_sample, S=SS, full = FALSE, N=200, no_short_sales = 1:14)
res_no_btc = EfficientFrontier(r= expected_return_sample[2:14], S=SS[2:14,2:14],full = FALSE, N =200, no_short_sales = 1:13)

lines(res_btc$sigma, res_btc$expected_return, col= 'green')
lines(res_no_btc$sigma[2:201], res_no_btc$expected_return[2:201], col = 'orange')

title(main = "Efficient Markowitz Mean Variance Frontier")
legend("bottomleft", legend = c("btc", "no btc", "btc no short sale","no btc no short sale"),
       col=c("darkgreen","red","green","orange"), lwd = 3, lty = c(1,1,1,1), cex=0.75)



# # target return
# target = 0.2
# 
# w = OptimalAllocation(r=expected_return_sample,S=SS, target_return = target)
# w_no_btc = OptimalAllocation(r=expected_return_sample[2:14],S=SS[2:14,2:14], target_return = target)
# 
# cbind(w, c(0,w_no_btc))


targets = seq(0.05,0.25, by = 0.01)
l=length(targets)
alloc_btc = zeros(14,l)
rownames(alloc_btc)=colnames(my_returns[,2*(1:14)])
alloc_no_btc = zeros(13,l)
rownames(alloc_no_btc)=colnames(my_returns[,2*(2:14)])

for(i in 1:l){
  alloc_btc[,i]= OptimalAllocation(r=expected_return_sample,S=SS, target_return = targets[i], no_short_sales = 1:14)
  alloc_no_btc[,i]= OptimalAllocation(r=expected_return_sample[2:14],S=SS[2:14,2:14], target_return = targets[i],no_short_sales = 1:13)
}

alloc_btc


# ggplot2 library
library(ggplot2)
library(RColorBrewer)

# plot(targets,alloc_btc[1,], type='l', col=rainbow(14)[1], ylim = c(-1,1))
# for(i in 2:14){
#   lines(targets, alloc_btc[i,], col= rainbow(14)[i])
# }


Asset = rep(rownames(alloc_btc), l)
Return = rep(targets, each = 14)
Values = drop(matrix(alloc_btc, nrow = 1))




colorRampPalette(brewer.pal(9, "Spectral"))(14)

data <- data.frame(Asset,Return,Values)
ggplot(data, aes(x=Return, y=Values, fill=Asset)) + 
  geom_area(alpha=1 , size=1, colour="black") +
  ggtitle("Markowitz Optimal Allocation Without Shortsellling")+
  scale_fill_manual(values =colorRampPalette(brewer.pal(9, "Paired"))(14) )





##################################################
with_btc = rbind(rep(0,length(alloc_no_btc)), alloc_no_btc)

with_btc[which(with_btc<1e-6)]=0

Asset = rep(rownames(alloc_btc), l)
Return = rep(targets, each = 14)
Values = drop(matrix(with_btc, nrow = 1))




colorRampPalette(brewer.pal(9, "Spectral"))(14)

data <- data.frame(Asset,Return,Values)
ggplot(data, aes(x=Return, y=Values, fill=Asset)) + 
  geom_area(alpha=1 , size=1, colour="black") +
  ggtitle("Markowitz Optimal Allocation Without Shortsellling")+
  scale_fill_manual(values =colorRampPalette(brewer.pal(9, "Paired"))(14) )

