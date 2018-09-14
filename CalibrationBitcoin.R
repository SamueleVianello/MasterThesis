
library(tseries)
library(mvtnorm)
library(pracma)
library(DEoptim)
library(statmod)
library(NMOF)


library(readxl)
library(xlsx)




library(binaryLogic)
library(tseries)


# Calibration on 3 asset including bitcoin

my_data = read.xlsx("XBT Correlations.xlsm", sheetName = "STATIC")



bitcoin = list(date = as.vector(my_data[1]), value =  as.vector(my_data[2]))
sp500 = list(date = as.vector(my_data[5]), value =  as.vector(my_data[6]))
  
  
plot(bitcoin$date[1:200,],bitcoin$value[1:200,],type='l')
plot(sp500$date[1:200,],sp500$value[1:200,],type = 'l',col='green')

  
btc_x = as.double(bitcoin$value[1:500,])
sp500_x = sp500$value[1:500,]

assets = cbind(log(btc_x[1:499]/btc_x[2:500]),log(sp500_x[1:499]/sp500_x[2:500]))
  
  
dt = 1/255

control_list = list(itermax = 2000, NP = 200, strategy = 6,trace=5)
### no common jump
bounds_nocommon = BoundsCreator(2, n_common=0)

start_time <- Sys.time()
outDE <- DEoptim(negloglik_2assets_nocommon,
                 lower = bounds_nocommon$lower,
                 upper = bounds_nocommon$upper,
                 control = control_list, dt = dt, x = assets, n=2)

end_time <- Sys.time()
end_time-start_time

ParametersReconstruction(outDE$optim$bestmem,2,common = FALSE)

