library(NMOF)
source("BatesModel.R")
source("CalibrationBates.R")



#########################################
####### SIMPLE PLOT TO SHOW DENSITIES ###
#########################################

r=-0.8
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

xx = seq(from=-10, to = 10, by =20/100)
all_dt = rep(t, length(xx))

yyH=pdfHeston(x=xx,x_0 = log(S_0),dt = all_dt,sigma_0 = sigma_0,r = r,k = k,eta = eta, theta = theta,rho = rho)
yyB=pdfBates(x=xx,x_0 = log(S_0),dt = all_dt,sigma_0 = sigma_0,r = r,k = k,eta = eta, theta = theta,rho = rho,
             lambda= lambda, mu_j = mu_j, sigma_j = sigma_j)

plot(xx,yyH,type='l',col='black')
grid()
lines(xx,yyB, col='blue')
sum(yyB*(xx[2]-xx[1]))
sum(yyH*(xx[2]-xx[1]))
legend("topright", legend = c("Heston","Bates"),col=c('black', 'blue'),
       lwd=3,lty=c(1,1),cex=0.75)


#########
# simple Test on loglikelihood

param = c(r,k,eta*0.8,22*theta,rho)

qx =rnorm(100, mean = 0.3)
qdt = rep(t,length(qx))
negloglikHeston(param, x=qx, x_0=log(S_0), sigma_0 = sigma_0, dt =qdt)

paramB=c(param,lambda,mu_j,sigma_j)
negloglikBates(paramB, x=qx, x_0=log(S_0), sigma_0 = sigma_0, dt =qdt)


#################################################################
## Test to see how eta and theta influence the distribution
#################################################################
n=20
theta_mix = seq(from=0.001, to=10, length.out = n)


xx = seq(from=-10, to = 10, by =20/100)
all_dt = rep(t, length(xx))

plot(0,0, ylim = c(0,0.7),xlim = c(-5,5))
grid()
for (i in 1:n){
  th =theta_mix[i]
  # print(th)
  yyH=pdfHeston(x=xx,x_0 = log(S_0),dt = all_dt,sigma_0 = sigma_0,r = r, k = k,eta = eta, theta = th,rho = rho)
  yyB=pdfBates(x=xx,x_0 = log(S_0),dt = all_dt,sigma_0 = sigma_0,r = r, k = k,eta = eta, theta = th,rho = rho,
               lambda =lambda, mu_j = mu_j, sigma_j = sigma_j)
  
  lines(xx,yyH)
  lines(xx,yyB,col='blue')
}

# the higher theta, the more leptokurtik it is




############

# calibration on gaussian generated data
set.seed(1234)
N_test=200
test_x = rnorm(N_test, mean = 0.01, sd = 0.6)

test_dt = rep(1, N_test)

test_sigma_0 = 0.3
test_x_0 = 0.02


initial = runif(5,min=-1)
params_H1 = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial, deoptim=FALSE)
#params$k*params$eta*2 > params$theta^2

initial = c(initial, runif(3,min=-1))
params_H2 = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial, deoptim=FALSE, model = "bates")
#paramsB = 

# hist( test_x, breaks = 20,freq = FALSE)
xx= seq(from=-3,to=3,by=1/50)

windows(width=10, height=8)
hist(test_x, freq = FALSE,breaks = seq(-3,3, length.out = 30))
lines(xx,dnorm(xx,mean=mean(test_x), sd=sd(test_x)), type='l', ylim = c(0,1))


yy1=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
             r = params_H1$r, k = params_H1$k, eta = params_H1$eta, theta = params_H1$theta, rho = params_H1$rho)
lines(xx,yy1,col='blue',type='l')
yy2=pdfBates(x=xx,x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
              r = params_H2$r, k = params_H2$k, eta = params_H2$eta, theta = params_H2$theta, rho = params_H2$rho,
             lambda = params_H2$lambda, mu_j = params_H2$mu_j, sigma_j = params_H2$sigma_j)
lines(xx,yy2,col='green')

legend("topright", legend = c("Gaussian sim","Heston","Bates"),col=c('black', 'blue','green'),
       lwd=3,lty=c(1,1,1),cex=0.75)
params_H1$objective_function
params_H2$objective_function

# n=100 0.1036277 0.6982232 0.0000100 0.0000100 0.7704672
# n=100 0.1036278 0.6982235 0.4014963 0.0000100 0.7704674

# weird curvy bates 0.1692278  1.6296841  0.1242986  0.6033132 -0.5870094  0.7570830  0.7912934  0.0000100
dev.off()

##################
# Why is bates CURVY??
xx= seq(from=-30,to=30,length.out = 1000)

yyb = pdfBates(x=xx, x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
               r = 0.1692278,  k = 1.6296841, eta =0.1242986, theta = 0.6033132,rho = -0.5870094,
               lambda = 0.7570830,mu_j = 0.7912934, sigma_j = 0.00001)

plot(xx,yyb, type = 'l', col='green')
grid()

sum(yyb*(xx[2]-xx[1]))




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






