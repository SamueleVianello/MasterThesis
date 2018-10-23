# volatility process setimation

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
N=length(btc)
dn =255*6

asset = eurostoxx

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
roll_sd = rolling_stdev(rev(asset)[1:dn],window=5,skip = 5)*sqrt(255)

lines((1:length(roll_sd))*5-4,roll_sd, col='blue',lwd=2)
#test = acf(roll_sd, lag.max = 1)

coef = arima(roll_sd, order = c(1,0,0))
#(coef)
a = coef$coef[1]
var_eps = coef$sigma2


dt= 5/255
beta= (a)/dt
delta = sqrt(var_eps/dt)

k= 2*(beta)
eta = delta**2 / k
theta = 2*delta

c(k=k,eta=eta,theta=theta)




# MONTHLY SPACED *** MONTHLY SPACED *** MONTHLY SPACED *** MONTHLY SPACED *** 
roll_sd = rolling_stdev(rev(asset)[1:dn],window=21,skip = 21)*sqrt(255)

plot((1:length(roll_sd))*5-4,roll_sd,type='l', col='blue',lwd=2)
#test = acf(roll_sd, lag.max = 1)

coef = arima(roll_sd, order = c(1,0,0))
#(coef)
a = coef$coef[1]
var_eps = coef$sigma2


dt= 5/255
beta= (a)/dt
delta = sqrt(var_eps/dt)

k= 2*(beta)
eta = delta**2 / k
theta = 2*delta

c(k=k,eta=eta,theta=theta)





# VASICEK 
roll_var = rolling_variance(rev(asset)[1:dn],window=5,skip = 5)*255


plot((1:length(roll_var))*5-4,roll_var,type='l', col='blue',lwd=2)
#test = acf(roll_sd, lag.max = 1)

coef_lm = lm(roll_var[-1]~roll_var[-length(roll_var)])
# summary(coef_lm)

b0 = coef_lm$coefficients[1]
b1 = coef_lm$coefficients[2]

dt= 5/255


k = (1-b1)/dt
eta = b0/k
theta = sqrt(var(coef_lm$residuals)/dt)


c(k=k,eta=eta,theta=theta/(sd(asset)*sqrt(255)))




