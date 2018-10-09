library(NMOF)
source("BatesModel.R")
source("CalibrationBates.R")
library(maxLik)



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





##########################################
# calibration on gaussian generated data
############################################
set.seed(1234)
N_test=1000
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
sum(yy1*(xx[2]-xx[1]))
sum(yy2*(xx[2]-xx[1]))   
dev.off()



##################
# Why is bates CURVY?? 

xx= seq(from=-10,to=10, length.out = 1000)

yyb = pdfBates(x=xx, x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
               r = 0.1692278,  k = 1.6296841, eta =0.1242986, theta = 0.6033132,rho = -0.5870094,
               lambda = 0.7570830, mu_j = 0.7912934, sigma_j = 0.00001)


# yyh = pdfHeston(x=xx, x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
#                r = 0.1692278,  k = 1.6296841, eta =0.1242986, theta = 0.6033132,rho = -0.5870094)
plot(xx,yyb, type = 'l', col='green')
grid()
# # Plot for different upper bounds in the integral
# for(i in (1:5)*10){
#   yyb = pdfBates(x=xx, x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
#                  r = 0.1692278,  k = 1.6296841, eta =0.1242986, theta = 0.6033132,rho = -0.5870094,
#                  lambda = 0.7570830, mu_j = 0.7912934, sigma_j = 0.00001, upper = i)
#   lines(xx,yyb,col=rainbow(10)[i%/%10])
# }
# lines(xx,yyh,col='blue')

sum(yyb*(xx[2]-xx[1]))

# increasing upper extremum does nothing to get rid of oscillations



xx = seq(from=-5,to=5, length.out = 10)

uu= seq(from=0,to=10, length.out = 1000)
ycfb=my_cfBates(uu,x_0 = test_x_0, v0 = test_sigma_0^2, tau = rep(test_dt[1],length(uu)),
                r = 0.1692278,  k = 1.6296841, vT =0.1242986, sigma = 0.6033132,rho = -0.5870094,
                lambda = 0.7570830, muJ = 0.7912934, vJ = 0.00001^2) * exp(-1i*xx[1]*uu)

plot(uu,Re(ycfb),type='l')
grid()
for(i in 2:length(xx)){
  ycfb=my_cfBates(uu,x_0 = test_x_0, v0 = test_sigma_0^2, tau = rep(test_dt[1],length(uu)),
                  r = 0.1692278,  k = 1.6296841, vT =0.1242986, sigma = 0.6033132,rho = -0.5870094,
                  lambda = 0.7570830, muJ = 0.7912934, vJ = 0.00001^2) * exp(-1i*xx[i]*uu)
  lines(uu,Re(ycfb),col=rainbow(10)[i])
  print(i)
}










param = c(r,k,eta,theta,rho,lambda, mu_j, sigma_j)

sigma_0 = sd(btc[1:25])*sqrt(255)
x_0 = log(my_data$BITCOIN[dim(my_data)[1]])

dn=10
negloglikHeston(params = param[1:5], x = cum_returns_btc[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn])

neglog=0
for(i in 1:10){
  neglog=neglog+negloglikHeston(params = param[1:5], x = cum_returns_btc[i], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[i])
}

neglog



# 














#**** GEN SA **** GEN SA **** GEN SA **** GEN SA **** GEN SA **** GEN SA **** GEN SA ****
#**** GEN SA **** GEN SA **** GEN SA **** GEN SA **** GEN SA **** GEN SA **** GEN SA ****
#**** GEN SA **** GEN SA **** GEN SA **** GEN SA **** GEN SA **** GEN SA **** GEN SA ****

library(GenSA)
initial = c(0.2, 0.9, sigma_0^2, 0.2, runif(1,min=-1))
dn=255*2
params_sa= GenSA(fn = negloglikHeston, lower = bounds$lower, upper = bounds$upper,par = initial,control = c(max.time=300, maxit=200,verbose=TRUE),
                  x = cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn])
params_sa$par

out_nlminb = nlminb(start = params_sa$par, objective = negloglikHeston, lower = bounds$lower, upper = bounds$upper,
                    x=cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0=sigma_0, dt = time_intervals[1:dn],
                    control=list(eval.max = 1000,iter.max = 100, trace = 10))
out_nlminb$par

plot(time_intervals[1:dn], cum_returns_eurostoxx[1:dn], type='l')


set.seed(1234)
N_test=1000
test_x = rnorm(N_test, mean = 0.01, sd = 0.6)

test_dt = rep(1, N_test)

test_sigma_0 = 0.3
test_x_0 = 0.02


initial = runif(5,min=-1)
params_H1 = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial, deoptim=TRUE, model = "heston_ab")

params_sa= GenSA(fn = negloglikHeston, lower = bounds$lower, upper = bounds$upper, par = initial,
                 control = c(max.time=180, maxit=20,verbose=TRUE,simple.function=TRUE),
                 x = test_x, x_0 = test_x_0, sigma_0 = test_sigma_0, dt = test_dt)
params_sa$par

out_nlminb = nlminb(start = params_sa$par, objective = negloglikHeston, lower = bounds$lower, upper = bounds$upper,
                    x=test_x, x_0 = test_x_0, sigma_0=test_sigma_0, dt = test_dt,
                    control=list(eval.max = 1000,iter.max = 100, trace = 1))

out_nlminb$par
par_2 = ParametersReconstruction(out_nlminb$par, model = "heston_ab")


xx= seq(from=-3,to=3,by=1/50)

windows(width=10, height=8)
hist(test_x, freq = FALSE,breaks = seq(-3,3, length.out = 30))
lines(xx,dnorm(xx,mean=mean(test_x), sd=sd(test_x)), type='l', ylim = c(0,1))

yy1=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
              r = params_H1$r, k = params_H1$k, eta = params_H1$eta, theta = params_H1$theta, rho = params_H1$rho)
lines(xx,yy1,col='blue',type='l')
yy2=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
             r = par_2$r, k = par_2$k, eta = par_2$eta, theta = par_2$theta, rho = par_2$rho)
lines(xx,yy2,col='green')


legend("topright", legend = c("Gaussian sim","Heston","Heston 2"),col=c('black', 'blue','green'),
       lwd=3,lty=c(1,1,1),cex=0.75)
params_H1$objective_function
#params_H2$objective_function



negloglikHeston(params = out_nlminb$par,
                  x = test_x, x_0 = test_x_0, sigma_0 = test_sigma_0, dt = test_dt)

# optimization using GenSA appears to be good but quite time consuming

# trying on real data:


#################################
###### Calibration on data ######
#################################



attach(my_returns)

# EUROSTOXX ** EUROSTOXX ** EUROSTOXX ** EUROSTOXX ** EUROSTOXX ** EUROSTOXX ** 
cum_returns_eurostoxx = cumsum(rev(eurostoxx))
time_intervals = rev(as.double((as.Date((btc_date)) - as.Date(btc_date[length(btc_date)]) + 1)/365 ))

N = length(cum_returns_eurostoxx)
dn=100


x_0 = log(my_data$EUROSTOXX50[N])
sigma_0 = sd( cum_returns_eurostoxx[1:25])
initial =  c(0.02, 0.4,0.3, 0.2, -0.3)

# choosing optimization 
params_1 = CalibrateModel(x = cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
                          trace = 1,model = "heston_ab", deoptim = TRUE, initial = initial)

bounds=BoundsCreator(n=1, model="heston_ab") 
params_SA = GenSA(fn = negloglikHeston, x = cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn], model="heston_ab",
                              par = initial,lower = bounds$lower, upper = bounds$upper,
                  control = c(max.time=180, maxit=20,verbose=TRUE,simple.function=TRUE))

params_SA$par
out_nlminb = nlminb(start = params_SA$par, objective = negloglikHeston, lower = bounds$lower, upper = bounds$upper,
                    x=cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0=sigma_0, dt =  time_intervals[1:dn],
                    control=list(eval.max = 1000,iter.max = 100, trace = 1))

out_nlminb$par
par_2 = ParametersReconstruction(out_nlminb$par, model = "heston_ab")
par_2

xx= seq(from=-3,to=3,by=1/50)
yy2=pdfHeston(x=xx,x_0 = x_0, sigma_0 = sigma_0, dt = rep( time_intervals[dn],length(xx)),
              r = par_2$r, k = par_2$k, eta = par_2$eta, theta = par_2$theta, rho = par_2$rho)
plot(xx,yy2,col='green',type='l')

plot(time_intervals[1:dn], cum_returns_eurostoxx[1:dn], type='l')




cum_returns_btc = cumsum(rev(btc))
time_intervals = rev(as.double((as.Date((btc_date)) - as.Date(btc_date[length(btc_date)]) + 1)/365 ))

sigma_0 = sd(btc[1:25])*sqrt(255)
x_0 = log(my_data$BITCOIN[dim(my_data)[1]])




N = length(cum_returns_btc)
dn=100

initial = c(0.2, 2, sigma_0^2, 0.2, runif(1,min=-1))




params = CalibrateModel(x = cum_returns_btc[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
                        trace = 1,model = "heston_ab", deoptim = TRUE, initial = initial)

params = fminsearch(f = negloglikHeston, x = cum_returns_btc[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],x0 = initial,maxiter = 100,tol = 1e-5)



plot(time_intervals[1:dn], cum_returns_btc[1:dn], type='l')
# 0.26159197 1.50004615 0.01424088 1.00000000 1.00000000











################################################
for(i in 1:100){
  if(i<=50) tag=TRUE
  else {tag = FALSE
      initial=runif(5)}
  params = CalibrateModel(x = cum_returns_btc[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
                          trace = 1,model = "heston", deoptim = tag, initial = initial)
  full_params[i,1]=params$r
  full_params[i,2]=params$k
  full_params[i,3]=params$eta
  full_params[i,4]=params$theta
  full_params[i,5]=params$rho
  full_params[i,6]=params$objective_function
}

colnames(full_params)=c('r','k','eta','theta','rho','loglik')


full_params_complete = matrix(rep(0,6*100), nrow = 100)

for(i in 1:100){
  if(i<=50) tag=TRUE
  else {tag = FALSE
  initial=runif(5)}
  params = CalibrateModel(x = cum_returns_btc, x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals,
                          trace = 1,model = "heston", deoptim = tag, initial = initial)
  full_params_complete[i,1]=params$r
  full_params_complete[i,2]=params$k
  full_params_complete[i,3]=params$eta
  full_params_complete[i,4]=params$theta
  full_params_complete[i,5]=params$rho
  full_params_complete[i,6]=params$objective_function
}
colnames(full_params_complete)=c('r','k','eta','theta','rho','loglik')




# new res1 0.07313141 4.84737910 0.00001000 0.98992618 1.00000000   823.3702
# new res2  5    5    2    1    1  911.002
# new res3 0.3477615 0.9736082 0.0000100 1.0000000 1.0000000    854.2386

plot(time_intervals, cum_returns_btc, type='l')


################## eurostoxx #########################
cum_returns_eurostoxx = cumsum(rev(eurostoxx))
time_intervals = rev(as.double((as.Date((btc_date)) - as.Date(btc_date[length(btc_date)]) + 1)/365 ))

full_params_eurostoxx = matrix(rep(0,6*100), nrow = 100)

# params = CalibrateModel(x = cum_returns_btc[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
#                         trace = 1,model = "heston", deoptim = TRUE, initial = initial)


for(i in 51:100){
  if(i<=50) tag=TRUE
  else {
    tag = FALSE
    initial=c(runif(4),runif(1,min=-1))
  }
  params = CalibrateModel(x = cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
                          trace = 1,model = "heston", deoptim = tag, initial = initial)
  full_params_eurostoxx[i,1]=params$r
  full_params_eurostoxx[i,2]=params$k
  full_params_eurostoxx[i,3]=params$eta
  full_params_eurostoxx[i,4]=params$theta
  full_params_eurostoxx[i,5]=params$rho
  full_params_eurostoxx[i,6]=params$objective_function
}

colnames(full_params_eurostoxx)=c('r','k','eta','theta','rho','loglik')


full_params_complete_eurostoxx = matrix(rep(0,6*100), nrow = 100)



params = CalibrateModel(x = cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
                        trace = 1,model = "heston", deoptim = TRUE, initial = initial)

for(i in 1:100){
  if(i<=50) tag=TRUE
  else {tag = FALSE
  initial=runif(5)}
  params = CalibrateModel(x = cum_returns_eurostoxx, x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals,
                          trace = 1,model = "heston", deoptim = tag, initial = initial)
  full_params_complete_eurostoxx[i,1]=params$r
  full_params_complete_eurostoxx[i,2]=params$k
  full_params_complete_eurostoxx[i,3]=params$eta
  full_params_complete_eurostoxx[i,4]=params$theta
  full_params_complete_eurostoxx[i,5]=params$rho
  full_params_complete_eurostoxx[i,6]=params$objective_function
}

colnames(full_params_complete_eurostoxx)=c('r','k','eta','theta','rho','loglik')

plot(time_intervals[1:dn], cum_returns_eurostoxx[1:dn], type='l')





