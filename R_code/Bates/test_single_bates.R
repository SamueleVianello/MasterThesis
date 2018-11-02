library(NMOF)
source("BatesModel.R")
source("CalibrationBates.R")
library(maxLik)



#########################################
####### SIMPLE PLOT TO SHOW DENSITIES ###
#########################################

r=0.1
t=1
S_0 = 1
sigma_0 = 0.05
theta =0.2
rho=0.9
k=2
eta = 1.5
lambda=5
mu_j = 0.1
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
negloglikHeston(param, x=qx, x_0=log(S_0), sigma_0 = sigma_0, dt =qdt, model = "heston", check_feller = FALSE)

paramB=c(param,lambda,mu_j,sigma_j)
negloglikBates(paramB, x=qx, x_0=log(S_0), sigma_0 = sigma_0, dt =qdt, model="bates")


#################################################################
## Test to see how eta and theta influence the distribution
#################################################################
n=10
theta_mix = seq(from=-1, to=1, length.out = n)


xx = seq(from=-10, to = 10, by =10/100)
all_dt = rep(t, length(xx))

plot(0,0, ylim = c(0,0.4),xlim = c(-5,5))
grid()
for (i in 1:n){
  th =(theta_mix)[i]
  # print(th)
  yyH=pdfHeston(x=xx,x_0 = log(S_0),dt = all_dt,sigma_0 = sigma_0,r = r, k = k,eta = eta, theta = theta,rho = th,unconditional = T)
  # yyB=pdfBates(x=xx,x_0 = log(S_0),dt = all_dt,sigma_0 = sigma_0,r = r, k = k,eta = eta,theta = th,rho = rho,
  #              lambda =lambda, mu_j = mu_j, sigma_j = sigma_j)
  
  print(paste("skewness =", 3/2*th*theta*sqrt(all_dt[1])/sqrt(eta)))
  lines(xx,yyH)
  #lines(xx,yyB,col='blue')
}

legend("topright", legend = c("Heston", "Bates"),col=c('black', 'blue'),
       lwd=2,lty=c(1,1),cex=0.75)

# the higher theta, the more leptokurtik it is


#################################################################
## Test to see how the pdf changes in time
#################################################################

r=1.5
t=1
S_0 = 1
sigma_0 = 0.1
theta =1.5
rho=0.9
k=2
eta = 0.4
lambda=5
mu_j = 0.1
sigma_j = 0.17


t= seq(from=0.1, to=1, by=0.1)


xx = seq(from=-10, to = 10, by =1/100)
all_dt = rep(t, length(xx))

plot(0,0, ylim = c(0,3),xlim = c(-5,5))
grid()
for (i in 1:length(t)){
  # print(th)
  yyH=pdfHeston(x=xx,x_0 = log(S_0),dt = rep(t[i], length(xx)), sigma_0 = sigma_0,r = r, k = k,eta = eta, theta = theta,rho = rho)
  yyB=pdfBates(x=xx,x_0 = log(S_0),dt = rep(t[i], length(xx)), sigma_0 = sigma_0,r = r, k = k,eta = eta, theta = theta,rho = rho,
               lambda =lambda, mu_j = mu_j, sigma_j = sigma_j)
  print(paste("Heston:", sum(yyH*(xx[2]-xx[1]))))
  print(paste("Bates:", sum(yyB*(xx[2]-xx[1]))))
  lines(xx,yyH)
  lines(xx,yyB,col='blue')
}
title(paste('sigma_0 =', sigma_0,'r =', r, 'k =', k,'eta =', eta, 'theta =', theta,'rho =', rho))
legend("topright", legend = c("Heston", "Bates"),col=c('black', 'blue'),
       lwd=2,lty=c(1,1),cex=0.75)


######################
### 3D plot in time

library(rgl)
# x = values to compute pdf
# y = time
# z = values of pdf


t= seq(from=0.05, to=1, by=0.05)

origin= matrix(c(0,0,0),ncol=3)
plot3d(x=origin, xlim = c(-5,5), ylim = c(min(t),max(t), zlim = c(0,2)),forceClipregion = TRUE)
aspect3d(c(1,1,0.1))

xx = seq(from=-5, to = 5, by =1/100)
for (i in 1:length(t)){
  # print(th)
  zzH=pdfHeston(x=xx,x_0 = log(S_0),dt = rep(t[i], length(xx)), sigma_0 = sigma_0,r = r, k = k,eta = eta, theta = theta,rho = rho)
  print(paste("Heston:", sum(zzH*(xx[2]-xx[1]))))
  
  yy= rep(t[i], length(xx))
  
  X = cbind(xx,yy,zzH)
  lines3d(X, zlim = c(0,2))
  # yyB=pdfBates(x=xx,x_0 = log(S_0),dt = rep(t[i], length(xx)), sigma_0 = sigma_0,r = r, k = k,eta = eta, theta = theta,rho = rho,
  #              lambda =lambda, mu_j = mu_j, sigma_j = sigma_j)
  # print(paste("Bates:", sum(yyB*(xx[2]-xx[1]))))
}




############################################
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












# GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN ***
# GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN ***
# GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN ***
# GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN *** GAUSSIAN ***
# test on gaussian generated values
set.seed(1234)
N_test=100
test_x = rnorm(N_test, mean = 0.01, sd = 0.6)

test_dt = rep(1, N_test)

test_sigma_0 = 0.3
test_x_0 = 0.02


initial = runif(5,min=-1)

start_time <- Sys.time()
params_H1 = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial, deoptim=TRUE, 
                           model = "heston", sigma_is_param = TRUE)
params_H2 = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial, deoptim=TRUE, 
                           model = "heston_ab", sigma_is_param = TRUE)
params_H1_noF = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial, deoptim=TRUE,
                               model = "heston", feller = TRUE, sigma_is_param = FALSE)
params_H2_noF = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial, deoptim=TRUE,
                               model = "heston_ab", feller = TRUE, sigma_is_param = FALSE)
end_time <- Sys.time()

end_time-start_time


xx= seq(from=-3,to=3,by=1/50)

windows(width=10, height=8)
hist(test_x, freq = FALSE,breaks = seq(-3,3, length.out = 30), ylim = c(0,1.2))
lines(xx,dnorm(xx,mean=mean(test_x), sd=sd(test_x)), type='l', ylim = c(0,1))

yy1=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = params_H1$sigma_0, dt = rep(test_dt[1],length(xx)),
              r = params_H1$r, k = params_H1$k, eta = params_H1$eta, theta = params_H1$theta, rho = params_H1$rho)
lines(xx,yy1,col='blue',type='l')
yy2=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = params_H2$sigma_0, dt = rep(test_dt[1],length(xx)),
             r = params_H2$r, k = params_H2$k, eta = params_H2$eta, theta = params_H2$theta, rho = params_H2$rho)
lines(xx,yy2,col='red')

yy3=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
              r = params_H1_noF$r, k = params_H1_noF$k, eta = params_H1_noF$eta, theta = params_H1_noF$theta, rho = params_H1_noF$rho)
lines(xx,yy3,col='aquamarine',type='l')

yy4=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
              r = params_H2_noF$r, k = params_H2_noF$k, eta = params_H2_noF$eta, theta = params_H2_noF$theta, rho = params_H2_noF$rho)
lines(xx,yy4,col='orange',type='l')


legend("topright", legend = c("Gaussian sim","Heston","Heston_ab","Heston_noF","Heston_ab_noF"),col=c('black', 'blue',"red","aquamarine", "orange"),
       lwd=3,lty=c(1,1,1),cex=0.75)
params_H1$objective_function
#params_H2$objective_function

res = cbind(c(params_H1$r,params_H1$k, params_H1$eta, params_H1$theta, params_H1$rho),
      c(params_H2$r,params_H2$k, params_H2$eta, params_H2$theta, params_H2$rho),
      c(params_H1_noF$r,params_H1_noF$k, params_H1_noF$eta, params_H1_noF$theta, params_H1_noF$rho),
      c(params_H2_noF$r,params_H2_noF$k, params_H2_noF$eta, params_H2_noF$theta, params_H2_noF$rho))
colnames(res)= c("Heston","Heston_ab","Heston_noF","Heston_ab_noF")
rownames(res)=c("r","k","eta","theta","rho")
res


negloglikHeston(c(params_H1$r,params_H1$k, params_H1$eta, params_H1$theta, params_H1$rho, params_H1$sigma_0),x = test_x,dt = test_dt,x_0=test_x_0, sigma_0 = test_sigma_0,
                model = "heston", check_feller = TRUE)



# test on a skewed distribution


set.seed(1234)
N_test=200
test_x = rexp(N_test,rate = 1.5)

test_dt = rep(0.5, N_test)

test_sigma_0 = 0.3
test_x_0 = 0.02


initial = runif(6,min=-1)

start_time <- Sys.time()
params_H1 = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial, deoptim=TRUE, 
                           model = "heston", sigma_is_param = TRUE)
params_H2 = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial, deoptim=TRUE, 
                           model = "heston_ab", sigma_is_param = TRUE)
params_H1_noF = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial[1:5], deoptim=TRUE,
                               model = "heston", feller = TRUE, sigma_is_param = FALSE)
params_H2_noF = CalibrateModel(x = test_x, x_0 = test_x_0,sigma_0 = test_sigma_0, dt = test_dt, trace = 1, initial[1:5], deoptim=TRUE,
                               model = "heston_ab", feller = TRUE, sigma_is_param = FALSE)
end_time <- Sys.time()

end_time-start_time


xx= seq(from=-3,to=5,by=1/50)

windows(width=10, height=8)
hist(test_x, freq = FALSE, breaks = seq(-2,5, length.out = 30), ylim = c(0,1.3))
lines(xx,dexp(xx,rate = 1.5), type='l')

yy1=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = params_H1$sigma_0, dt = rep(test_dt[1],length(xx)),
              r = params_H1$r, k = params_H1$k, eta = params_H1$eta, theta = params_H1$theta, rho = params_H1$rho)
lines(xx,yy1,col='blue',type='l')
yy2=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = params_H2$sigma_0, dt = rep(test_dt[1],length(xx)),
              r = params_H2$r, k = params_H2$k, eta = params_H2$eta, theta = params_H2$theta, rho = params_H2$rho)
lines(xx,yy2,col='red')

yy3=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
              r = params_H1_noF$r, k = params_H1_noF$k, eta = params_H1_noF$eta, theta = params_H1_noF$theta, rho = params_H1_noF$rho)
lines(xx,yy3,col='aquamarine',type='l')

yy4=pdfHeston(x=xx,x_0 = test_x_0, sigma_0 = test_sigma_0, dt = rep(test_dt[1],length(xx)),
              r = params_H2_noF$r, k = params_H2_noF$k, eta = params_H2_noF$eta, theta = params_H2_noF$theta, rho = params_H2_noF$rho)
lines(xx,yy4,col='orange',type='l')


legend("topright", legend = c("Exp sim","Heston","Heston_ab","Heston_no_sigma","Heston_ab_no_sigma"),col=c('black', 'blue',"red","aquamarine", "orange"),
       lwd=3,lty=c(1,1,1),cex=0.75)

sum(yy1*rep(xx[2]-xx[1], length(yy1)))
sum(yy2*rep(xx[2]-xx[1], length(yy1)))
sum(yy3*rep(xx[2]-xx[1], length(yy1)))
sum(yy4*rep(xx[2]-xx[1], length(yy1)))

res = cbind(c(params_H1$r,params_H1$k, params_H1$eta, params_H1$theta, params_H1$rho,params_H1$sigma_0, params_H1$objective_function),
            c(params_H2$r,params_H2$k, params_H2$eta, params_H2$theta, params_H2$rho,params_H2$sigma_0, params_H2$objective_function),
            c(params_H1_noF$r,params_H1_noF$k, params_H1_noF$eta, params_H1_noF$theta, params_H1_noF$rho,test_sigma_0, params_H1_noF$objective_function),
            c(params_H2_noF$r,params_H2_noF$k, params_H2_noF$eta, params_H2_noF$theta, params_H2_noF$rho,test_sigma_0, params_H2_noF$objective_function))
colnames(res)= c("Heston","Heston_ab","Heston_no_sigma","Heston_ab_no_sigma")
rownames(res)=c("r","k","eta","theta","rho","sigma_0","negloglik")
res


# trying on real data:


#################################
###### Calibration on data ######
#################################



attach(my_returns)


# EUROSTOXX ** EUROSTOXX ** EUROSTOXX ** EUROSTOXX ** EUROSTOXX ** EUROSTOXX ** 
cum_returns_eurostoxx = cumsum(rev(eurostoxx))
time_intervals = rev(as.double((as.Date((btc_date)) - as.Date(btc_date[length(btc_date)]) + 1)/365 ))

N = length(cum_returns_eurostoxx)
dn=255*6

# x_0 = log(my_data$EUROSTOXX50[N])
x_0 = 0
sigma_0 = sd( rev(eurostoxx)[1:dn])*sqrt(255)

# initial_mu = mean(cum_returns_eurostoxx[1:dn]/time_intervals[1:dn])
initial_mu = mean(rev(eurostoxx)[1:dn])*255
initial =  c(initial_mu, 4, sigma_0^2, 0.2, -0.3, sigma_0)


# res1=c(params_1$r,params_1$k, params_1$eta, params_1$theta, params_1$rho)
# res2=c(params_2$r,params_2$k, params_2$eta, params_2$theta, params_2$rho)
# for(i in 1:100){
t1=Sys.time()
params_1 = CalibrateModel(x = cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
                          trace = 1, model = "heston", deoptim = F, initial = initial, feller = F, sigma_is_param = TRUE)
# 0.21828180  0.87790103  0.03925139  0.26252176 -0.43957532  0.00001000 -219.2735
t2=Sys.time()
t2-t1

#                       alpha = k*eta       beta=k
initial =c(initial[1],initial[2]*initial[3],initial[2],initial[4],initial[5],initial[6]) 
t1=Sys.time()
params_2 = CalibrateModel(x = cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
                          trace = 1,model = "heston_ab", deoptim = F, initial = initial, feller= F, sigma_is_param = T)
t2=Sys.time()
t2-t1



# t1=Sys.time()
# params_B = CalibrateModel(x = cum_returns_eurostoxx[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
#                           trace = 1,model = "bates", deoptim = F, initial = c(initial,2,0.2,0.1) , feller= F, sigma_is_param = F)
# t2=Sys.time()
# t2-t1



# res1 = cbind(res1,c(params_1$r,params_1$k, params_1$eta, params_1$theta, params_1$rho, params_1$sigma_0, params_1$objective_function) )
# res2=cbind(res2,c(params_2$r,params_2$k, params_2$eta, params_2$theta, params_2$rho,sigma_0, params_2$objective_function))
# 
# }


# xx= seq(from=-30+x_0,to=30+x_0,by=1/10)
# windows(width=10, height=8)
# yy1=pdfHeston(x=xx,x_0 = x_0, sigma_0 = params_1$sigma_0, dt = rep( time_intervals[dn],length(xx)),
#               r = params_1$r, k = params_1$k, eta = params_1$eta, theta = params_1$theta, rho = params_1$rho)
# plot(xx,yy1,col='blue',type='l')
# 
# yy2 = pdfHeston(x=xx,x_0 = x_0, sigma_0 = sigma_0, dt = rep( time_intervals[dn],length(xx)),
#                 r = params_2$r, k = params_2$k, eta = params_2$eta, theta = params_2$theta, rho = params_2$rho)
# 
# lines(xx,yy2,col='green')
# 
# sum(yy1*rep(xx[2]-xx[1], length(yy1)))
# sum(yy2*rep(xx[2]-xx[1], length(yy1)))
# 
# res = cbind(c(params_1$r,params_1$k, params_1$eta, params_1$theta, params_1$rho,params_1$sigma_0, params_1$objective_function),
#             c(params_2$r,params_2$k, params_2$eta, params_2$theta, params_2$rho,sigma_0, params_2$objective_function))
# rownames(res)=c("r","k","eta","theta","rho","sigma_0","negloglik")
# res
# 
plot(time_intervals[1:dn], cum_returns_eurostoxx[1:dn], type='l')







cum_returns_sp500 = cumsum(rev(sp500))
time_intervals = rev(as.double((as.Date((sp500_date)) - as.Date(sp500_date[length(sp500_date)]) + 1)/365 ))


N = length(cum_returns_sp500)
dn=255*6

# x_0 = log(my_data$EUROSTOXX50[N])
x_0 = 0
sigma_0 = sd( rev(sp500)[1:dn])*sqrt(255)

# initial_mu = mean(cum_returns_eurostoxx[1:dn]/time_intervals[1:dn])
initial_mu = mean(rev(sp500)[1:dn])*255
initial =  c(initial_mu, 0.4, sigma_0^2, 0.2, -0.3, sigma_0)


t1=Sys.time()
params_1 = CalibrateModel(x = cum_returns_sp500[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = time_intervals[1:dn],
                          trace = 1, model = "heston", deoptim = F, initial = initial, feller = F, sigma_is_param = TRUE)
t2=Sys.time()
t2-t1


plot(time_intervals[1:dn], cum_returns_sp500[1:dn], type='l')






# Any given asset

asset = sp500
asset_date = eur_date

cum_returns_asset = cumsum(rev(asset))
time_intervals = rev(as.double((as.Date((asset_date)) - as.Date(asset_date[length(asset_date)]) + 1)/365 ))


N = length(cum_returns_asset)
dn=N


x_0 = 0
sigma_0 = sd( asset[1:dn])*sqrt(255)


initial =  c(mean(asset[1:dn])*255, 0.9, sigma_0^2, 0.01, -0.3, sigma_0)


initial[2]*initial[3]*2 > initial[4]^2


plot(time_intervals[1:dn], cum_returns_asset[1:dn], type='l')
set.seed(1234)

t1=Sys.time()
params_1 = CalibrateModel(x = asset[1:dn], x_0 = x_0, sigma_0 = sigma_0, dt = rep(1/255,dn),
                          trace = 1, model = "heston", deoptim = F, initial = initial, feller = T, sigma_is_param = F, unconditional = T)
t2=Sys.time()
t2-t1



params_1$k*params_1$eta*2 > params_1$theta^2

xx = seq(-.2,.2,length.out = 2000)
dt = rep(1/255,length(xx))

yy = pdfHeston(xx,x_0 = 0, dt = dt, sigma_0 = 100,r = params_1$r, k=params_1$k, eta = params_1$eta, theta = params_1$theta, rho = params_1$rho, unconditional = T)

hist(asset[1:dn],breaks = 100, freq = FALSE)
lines(xx,yy, type ='l', col='black')


# # add lines
# par = c(0.0672352, 0.675262, 0.0187067, 0.158945, -0.14)
# xx = seq(-0.4,0.4,length.out = 2000)
# dt = rep(1/255,length(xx))
# yy = pdfHeston(xx,x_0 = 0, dt = dt, sigma_0 = 100,r = par[1], k=par[2],eta = par[3],theta = par[4], rho = par[5], unconditional = T)
# hist(asset[1:dn],breaks = 100, freq = FALSE)
# lines(xx,yy, type ='l', col='blue')
# grid()


cdf_empiric = cumsum(yy*(xx[2]-xx[1]))

samples = interp1(cdf_empiric, xx, runif(100000, min = 0.0001, max=0.9999))
hist(samples,freq=F, breaks = 80)
lines(xx,yy, type ='l', col='blue')
skewness(samples)
