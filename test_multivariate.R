###################################################
############### Multivariate Calibration ##########
###################################################
rm(list = ls())

source("MultivariateMertonModel.R")


#####   chiamata delle librerie

library(tseries)
library(mvtnorm)
library(pracma)
library(DEoptim)
library(statmod)
library(NMOF)


library(pracma)
library(DEoptim)
library(statmod)
library(readxl)
library(xlsx)




library(binaryLogic)



##### testing the mvMertonpdf
#* (x, dt, mu, S, theta, delta, lambda, theta_z, delta_z, lambda_z, alpha)
n = 2
m = (c(0,0))
SS = matrix(c(0.8, 0.2,
              0.2, 0.1),ncol = 2)

thet = (c(0.5,1.1))
delt = (c(0.1,0.3))
lambd = (c(30,20))
thet_z = 1
delt_z = 0.5
lambd_z = 10
alph = (c(1,2))

xx = rbind(c(0,0))
dt = 1/255

dmvnorm(xx,mean = m, sigma = SS)

res =MultivariateMertonPdf(xx, dt=dt, m,SS,thet,delt,lambd,thet_z,delt_z, lambd_z, alph)
res



N=50
x <- seq(-5,5,length.out=N)
y <- seq(-5,5,length.out=N)

z <- matrix(rep(0,N*N),ncol = N)
X = matrix(rep(0,N*N),ncol = N)
Y = matrix(rep(0,N*N),ncol = N)
z_nojumps= matrix(rep(0,N*N),ncol = N)
for(i in 1:N){
  for(j in 1:N){
    X[i,j] = x[i]
    Y[i,j] = y[j]
    z[i,j] = MultivariateMertonPdf(c(x[i],y[j]), dt=dt, m,SS,thet,delt,lambd,thet_z,delt_z, lambd_z, alph)
    z_nojumps[i,j] = dmvnorm(c(x[i],y[j]),mean = m, sigma = SS)
  }
}



library(rgl)
surface3d(X,Y,z, col='green')
surface3d(X,Y, z_nojumps, col = 'blue')

mean((z-z_nojumps)^2)








################################ negloglikelihood function #############################

negloglik = function(params, x, dt, n) {
  # 
  # x is a matrix [Npoints * n] of all the points for which we compute the likelihood
  # 
  
  
  ## add check on inputs
  
  # reconstruction of parameters:
  idx =1
  mu = params[idx:(idx+n-1)]
  idx = idx+n 

  
  S = matrix(rep(0,n*n), ncol = n)
  i=1
  j=1
  for(k in 1:(n*(n+1)/2)){
    S[i,j] = params[idx+k-1]
    S[j,i] =  S[i,j]
    j=j+1
    if(j == n+1){
      i=i+1
      j=i
    }
  }
  idx = idx + n*(n+1)/2
  
  theta = params[idx:(idx+n-1)]
  idx = idx+n
  
  delta = params[idx:(idx+n-1)]
  idx = idx+n
  
  lambda = params[idx:(idx+n-1)]
  idx = idx+n
  
  theta_z = params[idx]
  idx = idx+1
  
  delta_z = params[idx]
  idx = idx+1
  
  lambda_z = params[idx]
  idx = idx+1
  
  alpha = params[idx:(idx+n-1)]
  idx = idx+n
  
  # print(mu)
  # print(S)
  # print(theta)
  # print(delta)
  # print(lambda)
  # print(alpha)
  if( (idx-1)!=length(params))
    stop("Error in parameter reconstruction: number of parameters is wrong.")
  
  
  # computing pdf on each point and adding
  partial = 0
  for(i in 1:dim(x)[1]){
    pdf = MultivariateMertonPdf(x[i,], dt, mu, S, theta, delta, lambda, theta_z, delta_z, lambda_z, alpha)
    # cat("\npdf:")
    # print(pdf)
    partial = partial + log(pdf)
  }
  
  # last check on results
  nll = -(partial)
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }
  return(nll)
}



## test negloglikelihood function

param=c(m,c(SS[1,1],SS[1,2],SS[2,2]), thet,delt,lambd,thet_z,delt_z,lambd_z,alph)
xx = rmvnorm(100, mean = m, sigma = SS)
negloglik(param, xx, dt=dt, n = 2)





###################### Calibration ####################

bounds = BoundsCreator(2, n_common=1)


outDE <- DEoptim(negloglik,
                 lower = bounds$lower,
                 upper = bounds$upper,
                 control = list(itermax = 100, NP = 100), dt = dt, x = xx, n=2)










## ==== DATA ====
## load data set (need a web connection)
x_EuroStoxx<-read.xlsx(file="EuroStoxx.xlsx",sheetName = "sheet1")[,2]
#x <- get.hist.quote(instrument = "AAPL",
# start = "2008-01-01", end = "2009-06-30",
# retclass = "zoo", quote = "AdjClose", compression = "d")



#x_EuroStoxx<-(simulateJump(0.1,0.25,1,-0.2,0.001,5,1/252)[,2])

#x_EuroStoxx<-ProcessoMerton(0.1,0.25,255*10,0.3,0.2,1/255)
## log-returns
dy <- diff(log(as.vector(x_EuroStoxx)))
## assume 255 days in a year (trading days)
dt <- 1 / 255
## ==== DEOPTIM ESTIMATION ====
#set.seed(1234)
outDE <- DEoptim(negloglik,
                 lower = c( 0.1, -10, 1e-4, -10, 1e-4),
                 upper = c(100, 10, 10, 10, 10),
                 control = list(itermax = 500, NP = 100), dt = dt, dy = dy)
#summary(outDE)
parametri=outDE$optim$bestmem

par=parametri[c(2,3,1,4,5)]
plot(x_EuroStoxx, type='l')
TotalTime=length(x_EuroStoxx)*dt

