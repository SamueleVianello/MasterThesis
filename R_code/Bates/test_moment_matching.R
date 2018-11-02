# test moment matching


library(moments)

source("CalibrationBates.R")

load("../returns.Rda")

attach(my_returns)

asset = sp500


skewness(asset)

kurtosis(asset)-3












mu = 0.0672352
k = 0.675262
eta = 0.0187067
theta =0.158945
rho = -0.141096
sigma_0 = 0.09414525
dt=1

2*k*eta > theta^2


xx = seq(-1,1,length.out = 1000)
dt = rep(1/255,length(xx))

yy = pdfHeston(xx,x_0 = 0, dt = dt, sigma_0 = sigma_0,r = mu, k=k,eta = eta,theta = theta,rho = rho, unconditional = T)

hist(sp500,breaks = 40, freq = FALSE)
lines(xx,yy, type ='l', col='blue')

