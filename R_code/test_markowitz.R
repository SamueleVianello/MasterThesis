source("MeanVariancePortfolio.R")

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


EfficientFrontier(r,S,full = FALSE)
