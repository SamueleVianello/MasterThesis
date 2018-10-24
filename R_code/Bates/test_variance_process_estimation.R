# volatility process estimation

rolling_stdev = function(data, window, skip=1){
  res= sd(data[1:window])
  
  idx=seq(1,length(data), skip)
  for(i in 2:(length(idx)-1)){
    res= c(res, sd(data[idx[i]:(idx[i]+window)]))
  }
  return(res)
}

rolling_variance = function(data, window, skip=1){
  res= var(data[1:window])
  
  idx=seq(1,length(data), skip)
  for(i in 2:(length(idx)-1)){
    res= c(res, var(data[idx[i]:(idx[i]+window)]))
  }
  return(res)
}


attach(my_returns)

asset = eurostoxx

N=length(asset)
dn =255*6



# DAILY SPACED *** DAILY SPACED *** DAILY SPACED *** DAILY SPACED ***
roll_sd = rolling_stdev(rev(asset)[1:dn],window=5,skip = 1)*sqrt(255)

plot(roll_sd, type='l', lwd=2)
#test = acf(roll_sd, lag.max = 1)

coef = arima(roll_sd, order = c(1,0,0))


# coef_lm = lm(roll_sd[-1]~roll_sd[-dn])
# summary(coef_lm)

a = coef$coef[1]
var_eps = coef$sigma2

dt= 1/255
beta= (1-a)/dt
delta = sqrt(var_eps/dt)

k= 2*(beta)
eta = delta**2 / k
theta = 2*delta

c(k=k,eta=eta,theta=theta)



# WEEKLY SPACED *** WEEKLY SPACED *** WEEKLY SPACED *** WEEKLY SPACED ***
days = 5
roll_sd = rolling_stdev(rev(asset)[1:dn],window=days,skip = days)*sqrt(255)

lines((1:length(roll_sd))*days-days+1,roll_sd, col='blue',lwd=2)
#test = acf(roll_sd, lag.max = 1)

coef = arima(roll_sd, order = c(1,0,0))
#(coef)
a = coef$coef[1]
var_eps = coef$sigma2


dt= days/255
beta= (1-a)/dt
delta = sqrt(var_eps/dt)

k= 2*(beta)
eta = delta**2 / k
theta = 2*delta

c(k=k,eta=eta,theta=theta)




# MONTHLY SPACED *** MONTHLY SPACED *** MONTHLY SPACED *** MONTHLY SPACED *** 
days = 21
roll_sd = rolling_stdev(rev(asset)[1:dn],window=days,skip = days)*sqrt(255)

lines((1:length(roll_sd))*days-days+1,roll_sd,type='l', col='green',lwd=2)
#test = acf(roll_sd, lag.max = 1)

coef = arima(roll_sd, order = c(1,0,0))
#(coef)
a = coef$coef[1]
var_eps = coef$sigma2


dt= days/255
beta= (1-a)/dt
delta = sqrt(var_eps/dt)

k= 2*(beta)
eta = delta**2 / k
theta = 2*delta

c(k=k,eta=eta,theta=theta)





# VASICEK 
days = 5
roll_var = rolling_variance(rev(asset)[1:dn],window=days,skip = 1)*255


plot((1:length(roll_var))*days-days+1,roll_var,type='l', col='blue',lwd=2)
#test = acf(roll_sd, lag.max = 1)

y=roll_var[-1]
x=roll_var[-length(roll_var)]
coef_lm = lm(y~x)
# summary(coef_lm)

points((1:length(coef_lm$fitted.values))*days-days+1,coef_lm$fitted.values, pch = 3)

b0 = coef_lm$coefficients[1]
b1 = coef_lm$coefficients[2]

dt= 1/255


k = (1-b1)/dt
eta = b0/(k*dt)
theta = sqrt(var(coef_lm$residuals)/dt)

c(b0 = b0, b1=b1, var_res = var(coef_lm$residuals), k=k,eta=eta,theta=theta/(sd(asset)*sqrt(255)))


