# Compare fft pdf with other pdf inversion

source("BatesModel.R")
source("AuxiliaryFunctions.R")
load("returns.Rda")

attach(my_returns)

#  some parameters calibrated on sp500 
mu = 0.127668393
k = 1.98833
eta = 0.020762953
theta =0.28734509
rho = -0.000107735

dt=1

check_feller(c(mu,k,eta,theta,rho))


xx = seq(-0.2,0.2,length.out = 10000)
dt = rep(1/255,length(xx))

yy = pdfHeston(xx, dt = dt,r = mu, k=k,eta = eta,theta = theta,rho = rho)

hist(sp500,breaks = 50, freq = FALSE)
points(xx,yy, type ='l', col='blue')





yy_fft_12 = pdfHeston_fft(x=xx, dt=1/255, r = mu,k = k,eta = eta,theta = theta,rho = rho,N = 2**12)

plot(xx, yy, type = 'l')
points(xx, yy_fft_12, col='blue',type='l',lwd=1)
grid()
legend('topright', legend=c('FFT', 'quadrature'), col =c('blue','black'), lty = 1 )


########################## bates ############################

mu = 0.127668393
k = 1.55897463
eta = 0.017922414
theta =0.236392
rho = -0.000347238

muj = -0.04555365		
sigmaj = 0.00001
lambdaj= 0.6237174


xx = seq(-0.7,0.7,length.out = 5000)
dt = rep(1/255,length(xx))

yy = pdfBates(xx, dt = dt,r = mu, k=k,eta = eta,theta = theta,rho = rho,
                     mu_j = muj,sigma_j = sigmaj, lambda = lambdaj)

hist(sp500,breaks = 100, freq = FALSE)
points(xx,yy, type ='l', col='blue')


yy_fft_bates = pdfBates(xx, dt = 1/255,r = mu, k=k,eta = eta,theta = theta,rho = rho,
                        mu_j = muj,sigma_j = sigmaj, lambda = lambdaj, N = 2**12)


plot(xx,yy_fft_bates, type='l')
lines(xx,yy, type ='l', col='blue')
grid()
legend('topright', legend=c('FFT', 'quadrature'), col =c('blue','black'), lty = 1 )

