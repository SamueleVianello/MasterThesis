source("CalibrationBates.R")
source("BatesModel.R")
source('UtilityFunctions.R')
library(moments)


# pdf bates

r = 0.12
k = 1
eta = 0.02
theta =0.1
rho = 0

mu_j = 0.1
sigma_j = 0.05
lambda = 1

sigma_0 = 0.12

check_feller(c(r,k,eta,theta,rho))

xx = seq(-1,1, length.out = 200)
dt = 1


yy_B = pdfBates(x=xx, dt = dt, r = r,k = k,eta = eta,theta = theta,rho = rho,
                       lambda = lambda, mu_j = mu_j, sigma_j = sigma_j)
yy_H = pdfHeston(x=xx, dt=dt, r = r, k = k,eta = eta,theta = theta,rho = rho)



plot(xx,yy_H, pch=3)
lines(xx, yy_B0, col='blue')
lines(xx, yy_B, col = 'red')
grid()


#### how parameters affect shapes #######


# changing lambda
l = (1:10)*0.5
cols = terrain.colors(length(l))

xx = seq(-1,1, length.out = 200)
dt = rep(1, length(xx))
yy_B = pdfBates(x=xx, dt = dt, r = r,k = k,eta = eta,theta = theta,rho = rho,
                       lambda = 0, mu_j = mu_j, sigma_j = sigma_j)

plot(xx,yy_B, type='l', main= 'Lambda')
for(i in 1:length(l)){
  yy=pdfBates(x=xx, dt = dt, r = r,k = k,eta = eta,theta = theta,rho = rho,
                       lambda = l[i], mu_j = mu_j, sigma_j = sigma_j)
  print(c(lambda = l[i],mean = mean_from_pdf(xx,yy),var = variance_from_pdf(xx,yy), skew =skewness_from_pdf(xx,yy), kurtosis = kurtosis_from_pdf(xx,yy)))
  lines(xx,yy, col= cols[i])
}


# changing mu_j
m_j = seq(-0.20, 0.2, length.out = 21)
cols = topo.colors(length(m_j))

xx = seq(-4,4, length.out = 500)
dt = rep(5, length(xx))
yy_B = pdfBates(x=xx, dt = dt, r = r,k = k,eta = eta,theta = theta,rho = rho,
                       lambda = lambda, mu_j = m_j[1], sigma_j = sigma_j)

plot(xx,yy_B, type='l',main= 'mu_j', ylim = c(0, 1.2))
for(i in 1:length(m_j)){
  yy=pdfBates(x=xx, dt = dt, r = r,k = k,eta = eta,theta = theta,rho = rho,
                     lambda = lambda, mu_j = m_j[i], sigma_j = sigma_j,
                     upper = 1000)
  print(c(mu_j= m_j[i],mean = mean_from_pdf(xx,yy),var = variance_from_pdf(xx,yy),skew =skewness_from_pdf(xx,yy), kurtosis = kurtosis_from_pdf(xx,yy),int_pdf = sum(yy*(xx[2]-xx[1]))))
  lines(xx,yy, col= cols[i], type='l')
}


# changing sigma_j
sig_j = seq(0.0, 1, length.out = 21)
cols = topo.colors(length(sig_j))

xx = seq(-4,4, length.out = 500)
dt = rep(5, length(xx))
yy_B = pdfBates(x=xx, dt = dt, r = r,k = k,eta = eta,theta = theta,rho = rho,
                       lambda = lambda, mu_j = mu_j, sigma_j = sig_j[1])

plot(xx,yy_B, type='l',main= 'mu_j', ylim = c(0, 1.2))
for(i in 1:length(sig_j)){
  yy=pdfBates(x=xx, dt = dt, r = r,k = k,eta = eta,theta = theta,rho = rho,
                     lambda = lambda, mu_j = mu_j, sigma_j = sig_j[i],
                     upper = 1000)
  print(c(sigma_j= sig_j[i],mean = mean_from_pdf(xx,yy),var = variance_from_pdf(xx,yy),skew =skewness_from_pdf(xx,yy), kurtosis = kurtosis_from_pdf(xx,yy),int_pdf = sum(yy*(xx[2]-xx[1]))))
  lines(xx,yy, col= cols[i], type='l')
}



