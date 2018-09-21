source("PPgen.R")
source("callMerton.R")
source("Allocazione.R")
source("ProcessoMerton.R")
source("ProcessoMertonCorrelazione.R")
source("PortfolioInsurance.R")


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

#############################################################################
###########################################################################
###### Calibrazione Jump Diffusion
##########################################################################
##############################################################################


## ==== MIXTURE DENSITY AND LIKELIHOOD ====
## mixture density
fdy <- function(dy, dt, lambda, mu, sigma, muq, sigmaq) {
  mu1 <- (mu - sigma^2 / 2) * dt
  mu2 <- (mu - sigma^2 / 2) * dt + muq
  sig1 <- sigma * sqrt(dt)
  sig2 <- sqrt(sigma^2 * dt + sigmaq^2)
  pdf1 <- dnorm(dy, mean = mu1, sd = sig1)
  pdf2 <- dnorm(dy, mean = mu2, sd = sig2)
  pdf <- (1 - lambda * dt) * pdf1 + (lambda * dt) * pdf2
  return(pdf)
}

## negloglikelihood function
negloglik <- function(theta, dy, dt) {
  L <- fdy(dy = dy, dt = dt, lambda = theta[1],
           mu = theta[2], sigma = theta[3],
           muq = theta[4], sigmaq = theta[5])
  nll <- -sum(log(L))
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll <- 1e10
  }
  return(nll)
}


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

# matplot(simulateJump(par[1],par[2],par[3],par[4],par[5],TotalTime,dt)[,2],type='l')
# matplot(simulateJump(par[1],par[2],par[3],par[4],par[5],TotalTime,dt)[,2],type='l')
# matplot(simulateJump(par[1],par[2],par[3],par[4],par[5],TotalTime,dt)[,2],type='l')
# matplot(simulateJump(par[1],par[2],par[3],par[4],par[5],TotalTime,dt)[,2],type='l')

# parametri EurostoXX

drift_1 = par[1]
sigma_1= par[2]
lambda_1=par[3]
mu_q_1=par[4]
sigma_q_1=par[5]
hist(dy,breaks=200,prob=TRUE,main="EuroStoxx",xlab="daily returns")
curve(dnorm(x, mean=mean(dy), sd=sd(dy)), add=TRUE,col="blue")
lines(sort(dy),fdy (dy=sort(dy),dt=dt, lambda=lambda_1, mu=drift_1, sigma=sigma_1,muq= mu_q_1, sigmaq=sigma_q_1),col="red")
#plot(fdy (dy, dt, lambda_1, drift_1, sigma_1, mu_q_1, sigma_q_1))
#
legend("topright", legend = c("Real", "Black-Scholes","Merton"),col=c('black','blue','red'),lwd=3, pch = c(NA, NA,NA),lty=c(1,1,1))
# matplot(sort(dy),fdy (dy=sort(dy),dt=dt, lambda=lambda_1, mu=drift_1, sigma=sigma_1,muq= mu_q_1, sigmaq=sigma_q_1),type='l', xlab="daily return",ylab="density")