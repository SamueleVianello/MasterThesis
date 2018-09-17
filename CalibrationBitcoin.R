
library(tseries)
library(mvtnorm)
library(pracma)
library(DEoptim)
library(statmod)
library(NMOF)


library(readxl)
library(xlsx)




library(binaryLogic)

source('MultivariateMertonModel.R')

# Calibration on 3 asset including bitcoin

my_data<-read.xlsx(file="XBT Correlations.xlsm",sheetName = "STATIC" , header=TRUE)


bitcoin = list(date = (my_data[[1]]), value =  (my_data[[2]]))
sp500 = list(date = (my_data[[5]]), value = (my_data[[6]]))

  
N=500

btc_x =(bitcoin$value[1:N])
sp500_x = sp500$value[1:N]
eurostoxx_euro_x =my_data$EUROSTOXX50[1:N] 


assets_return = cbind(log(btc_x[1:(N-1)]/btc_x[2:N]),
                      log(sp500_x[1:(N-1)]/sp500_x[2:N]),
                      log(eurostoxx_euro_x[1:(N-1)]/eurostoxx_euro_x[2:N]))

x11()
par(mfrow = c(2,3))
plot(bitcoin$date[1:N],bitcoin$value[1:N],type='l')
plot(sp500$date[1:N],sp500$value[1:N],type = 'l',col='green')
plot(my_data$EUROSTOXX50_DATE[1:N],my_data$EUROSTOXX50[1:N],type = 'l',col='blue')
plot(bitcoin$date[1:(N-1)],assets_return[,1], type='l')
plot(sp500$date[1:(N-1)],assets_return[,2], type='l')
plot(my_data$EUROSTOXX50_DATE[1:(N-1)],assets_return[,3], type='l',col='blue')

graphics.off()
  










############# CALIBRATION ######################

dt = 1/255

control_list = list(itermax = 3000, NP = 200, strategy = 6,trace=5)
### no common jump
bounds_nocommon = BoundsCreator(3, n_common=0)

start_time <- Sys.time()
outDE <- DEoptim(negloglik_3assets_nocommon,
                 lower = bounds_nocommon$lower,
                 upper = bounds_nocommon$upper,
                 control = control_list, dt = dt, x = assets_return, n=3)

end_time <- Sys.time()
calibration_time = end_time-start_time
calibration_time

calibrated = ParametersReconstruction(outDE$optim$bestmem,3,common = FALSE)
calibrated

correlation = cov2cor(calibrated$S)
correlation


#
negloglik_nocommon(outDE$optim$bestmem,assets_return,dt,n=3)
MultivariateMertonPdf_nocommon(assets_return[1,],dt=dt, mu = calibrated$mu, S = calibrated$S, 
                               theta = calibrated$theta, delta = calibrated$delta, lambda = calibrated$lambda)

