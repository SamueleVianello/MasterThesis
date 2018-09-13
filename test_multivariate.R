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
x <- seq(-1,1,length.out=N)
y <- seq(-1,1,length.out=N)

z <- matrix(rep(0,N*N),ncol = N)
X = matrix(rep(0,N*N),ncol = N)
Y = matrix(rep(0,N*N),ncol = N)
z_nojumps= matrix(rep(0,N*N),ncol = N)
for(i in 1:N){
  for(j in 1:N){
    X[i,j] = x[i]
    Y[i,j] = y[j]
    z[i,j] = MultivariateMertonPdf(c(x[i],y[j]), dt=dt, m,SS,thet,delt,lambd,thet_z,delt_z, lambd_z, alph)
    z_nojumps[i,j] = dmvnorm(c(x[i],y[j]),mean = m*dt, sigma = SS*sqrt(dt))
  }
}



library(rgl)
plot3d(X,Y,z, col='green')
points3d(X,Y, z_nojumps, col = 'blue')
legend3d("bottomleft",c("Jumps", "No jumps"))

mean((z-z_nojumps)^2)








################################ negloglikelihood function #############################

## test negloglikelihood function

param=c(m,c(SS[1,1],SS[1,2],SS[2,2]), thet,delt,lambd,thet_z,delt_z,lambd_z,alph)
xx = rmvnorm(1000, mean = m*dt, sigma = SS*dt)

start_time <- Sys.time()
negloglik(param, xx, dt=dt, n = 2)
end_time <- Sys.time()

end_time-start_time

lx =lapply(seq_len(nrow(xx)), function(i) xx[i,])

start_time <- Sys.time()
vnegloglik(param, lx, dt=dt, n = 2)
end_time <- Sys.time()

end_time-start_time

start_time <- Sys.time()
negloglik_2assets(param, xx, dt=dt, n = 2)
end_time <- Sys.time()

end_time-start_time



#######################################################
###################### Calibration ####################
control_list = list(itermax = 1000, NP = 200, strategy = 6,trace=5)

bounds = BoundsCreator(2, n_common=1)

start_time <- Sys.time()
outDE <- DEoptim(negloglik_2assets,
                 lower = bounds$lower,
                 upper = bounds$upper,
                 control = control_list, dt = dt, x = xx, n=2)

end_time <- Sys.time()
end_time-start_time


# 
# cov_z = alph%*%t(alph)*delt_z
# mean_z = thet_z*alph
# 
# dmvnorm(x, mean = m*dt + mean_z, sigma = SS*dt + cov_z)





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

