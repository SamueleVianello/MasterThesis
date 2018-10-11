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

expected_return = 0.04
  
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


EfficientFrontier(r,S,full = TRUE)


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
#### test on our  calibrated assets ######
##########################################
attach(my_returns)

# expected_return_sample = colMeans(my_returns[,2*(1:14)])

# yearly
expected_return_sample =  colSums(results$full_mu)/6 + colSums(results$full_theta)/6 * colSums(results$full_lambda)/6 

colSums(results$full_mu)/6 * 


#SS = cov(my_returns[,2*(1:14)])
SS = results$covariance

res1 = EfficientFrontier(expected_return_sample,SS,full = FALSE, plot = FALSE)
res2 = EfficientFrontier(expected_return_sample[2:14],SS[2:14,2:14],full = FALSE, plot = FALSE)

plot(res1$sigma,res1$expected_return, type = 'l',col='green', ylim = c(min(c(res1$expected_return,expected_return_sample)), max(res1$expected_return)))
lines(res2$sigma,res2$expected_return, col = 'red')
points(sqrt(diag(SS)),expected_return_sample, pch='+', col = 'blue')
text(sqrt(diag(SS)),expected_return_sample, labels = colnames(my_returns[,2*(1:14)]),pos = 3)
grid()

target = 0.002

w = OptimalAllocation(expected_return_sample,SS, sd = 0.04)
w_no_btc = OptimalAllocation(expected_return_sample[2:14],SS[2:14,2:14], sd = 0.04)

cbind(w, c(0,w_no_btc))

