library(NMOF)
source("BatesModel.R")
source("CalibrationBates.R")
u = 1

r=0.8
t=1
S_0 = 1
sigma_0 = 0.05
theta =0.8
rho=0.9
k=2
eta = 1.5
lambda=5
mu_j = -0.1
sigma_j = 0.17


# FELLER CONDITION
# 2*k*eta > theta^2




xx = seq(from=-10, to = 10, by =20/1000)
all_dt = rep(t, length(xx))


yyH=pdfHeston(x=xx,x_0 = log(S_0),dt = all_dt,sigma_0 = sigma_0,r = r,k = k,eta = eta,theta = theta,rho = rho)
yyB=pdfBates(x=xx,x_0 = log(S_0),dt = all_dt,sigma_0 = sigma_0,r = r,k = k,eta = eta,theta = theta,rho = rho, 
             lambda = lambda, mu_j = mu_j, sigma_j = sigma_j)

plot(xx,yyH,type='l')
grid()
lines(xx,yyB,col='blue')
sum(yyB*(xx[2]-xx[1]))
sum(yyH*(xx[2]-xx[1]))
#legend("topright", legend=())

#######

param = c(r,k,eta,theta,rho)

qx =rnorm(1000, mean = 0.3)
negloglikHeston(param, x=qx, x_0=log(S_0), sigma_0 = sigma_0, dt =all_dt[1:length(qx)])
paramB=c(param,lambda,mu_j,sigma_j)
negloglikBates(paramB, x=qx, x_0=log(S_0), sigma_0 = sigma_0, dt =all_dt[1:length(qx)])




############

# calibration on gaussian generated data
N_test=30
test_x = rnorm(N_test, mean = 0.01, sd = 0.6)

test_dt = rep(1, N_test)

test_sigma_0 = 0.3
test_x_0 = 0.02

set.seed(1234)
params = CalibrateHeston(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1)
params$k*params$eta*2 > params$theta^2

paramsB = 

# hist( test_x, breaks = 20,freq = FALSE)
plot(xx,dnorm(xx,mean=0.01,sd=0.6), type='l')
yy=pdfHeston(x=test_x,x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(test_x)),
             r = params$r, k = params$k, eta = params$eta, theta = params$theta, rho = params$rho)
plot(test_x,yy)







# 
# attach(my_returns)
# cum_returns_btc = cumsum(btc)
# time_intervals = as.double((as.Date(btc_date) - as.Date(btc_date[length(btc_date)]) + 1)/365 )
# 
# sigma_0 = sd(btc[1:25])*sqrt(255)
# x_0 = log(my_data$BITCOIN[dim(my_data)[1]])
# 
# param = c(r,k,eta,theta,rho)
# N = length(cum_returns_btc)
# dn=10
# negloglikHeston(param, x=cum_returns_btc[(N-dn):N],x_0=x_0, sigma_0 = sigma_0, dt = time_intervals[(N-dn):N])
# 
# 
# 
# params = CalibrateHeston(x = cum_returns_btc[(N-dn):N], x_0 = x_0,sigma_0 = sigma_0,dt = time_intervals[(N-dn):N],trace = 1)  
# 
# 



#COME STIMARE VARIANZA??






