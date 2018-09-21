

library(mvtnorm)
library(DEoptim)
library(statmod)
library(NMOF)


library(readxl)
library(xlsx)



library(binaryLogic)

source('MultivariateMertonModel.R')
source('CalibrationMVMerton.R')


# CREATION OF THE DATASET OF LOG-RETURNS FROM EXCEL DATASET

# my_data<-read.xlsx(file="XBT Correlations.xlsm",sheetName = "STATIC" , header=TRUE)
# 
# 
# 
# leng = dim(my_data)[1]
# my_returns = data.frame(btc_date = my_data$BTC_DATE[1:(leng-1)])
# my_returns$btc = log(my_data$BITCOIN[1:(leng-1)]/my_data$BITCOIN[2:leng])
# 
# my_returns$bric_date = my_data$MSCI.BRIC._DATE[1:(leng-1)]
# my_returns$bric = log(my_data$MSCI.BRIC.[1:(leng-1)]/my_data$MSCI.BRIC.[2:leng])
# 
# my_returns$sp500_date = my_data$S.P500_DATE[1:(leng-1)]
# my_returns$sp500 = log(my_data$S.P500[1:(leng-1)]/my_data$S.P500[2:leng])
# 
# my_returns$eurostoxx_date = my_data$EUROSTOXX50_DATE[1:(leng-1)]
# my_returns$eurostoxx = log(my_data$EUROSTOXX50[1:(leng-1)]/my_data$EUROSTOXX50[2:leng])
# 
# my_returns$gold_date = my_data$GOLD_DATE[1:(leng-1)]
# my_returns$gold = log(my_data$GOLD[1:(leng-1)]/my_data$GOLD[2:leng])
# 
# my_returns$wti_date = my_data$WTI_DATE[1:(leng-1)]
# my_returns$wti = log(my_data$WTI[1:(leng-1)]/my_data$WTI[2:leng])
# 
# my_returns$grain_date = my_data$GRAIN_DATE[1:(leng-1)]
# my_returns$grain = log(my_data$GRAIN[1:(leng-1)]/my_data$GRAIN[2:leng])
# 
# my_returns$metal_date = my_data$IND.METALS_DATE[1:(leng-1)]
# my_returns$metal = log(my_data$IND.METALS[1:(leng-1)]/my_data$IND.METALS[2:leng])
# 
# my_returns$eur_date = my_data$EUR_DATE[1:(leng-1)]
# my_returns$eur = log(my_data$EUR[1:(leng-1)]/my_data$EUR[2:leng])
# 
# my_returns$gbp_date = my_data$GBP_DATE[1:(leng-1)]
# my_returns$gbp = log(my_data$GBP[1:(leng-1)]/my_data$GBP[2:leng])
# 
# my_returns$chf_date = my_data$CHF_DATE[1:(leng-1)]
# my_returns$chf = log(my_data$CHF[1:(leng-1)]/my_data$CHF[2:leng])
# 
# my_returns$jpy_date = my_data$JPY_DATE[1:(leng-1)]
# my_returns$jpy = log(my_data$JPY[1:(leng-1)]/my_data$JPY[2:leng])
# 
# my_returns$pan_euro_date = my_data[[25]][1:(leng-1)]
# my_returns$pan_euro = log(my_data$BBG.Barclays.PAN.EURO.Aggregate[1:(leng-1)]/my_data$BBG.Barclays.PAN.EURO.Aggregate[2:leng])
# 
# my_returns$pan_us_date = my_data[[27]][1:(leng-1)]
# my_returns$pan_us = log(my_data$BBG.Barclays.PAN.US.Aggregate[1:(leng-1)]/my_data$BBG.Barclays.PAN.US.Aggregate[2:leng])
# 
#
#
# save(my_returns, file = "returns.Rda")

load("returns.Rda")
load("data.Rda")


N=500

# Plot of first few data
x11()
par(mfrow = c(2,3))
plot(my_data$BTC_DATE[1:N],my_data$BITCOIN[1:N],type='l')
plot(my_data$S.P500_DATE[1:N],my_data$S.P500[1:N],type = 'l',col='green')
plot(my_data$EUROSTOXX50_DATE[1:N],my_data$EUROSTOXX50[1:N],type = 'l',col='blue')
plot(my_returns$btc_date[1:N],my_returns$btc[1:N], type='l')
plot(my_returns$sp500_date[1:N],my_returns$sp500[1:N], type='l',col='green')
plot(my_returns$eurostoxx_date[1:N],my_returns$eurostoxx[1:N], type='l',col='blue')

graphics.off()
  




################ CALIBRATION USING deoptim + nlminb ############################
attach(my_returns)


N = dim(my_returns)[1]

dt = 1/255


# calibration_time = end_time-start_time
# calibration_time
# 
# calibrated = ParametersReconstruction(outDE$optim$bestmem,4,common = FALSE)
# calibrated
# 
# correlation = cov2cor(calibrated$S)
# correlation
# 
# 
# initial=outDE$optim$bestmem
# 
# bounds_nocommon = BoundsCreator(4, n_common=0)
# 
#
# start_time <- Sys.time()
# out_nlminb = nlminb(start = c(calibrated_params),objective = negloglik_4assets_nocommon,lower = bounds_nocommon$lower,
#        upper = bounds_nocommon$upper,dt=dt, x=cbind(btc[1:N],sp500[1:N],bric[1:N],eurostoxx[1:N]),n=4,
#        control=list(eval.max = 10000,iter.max = 1000, trace = 10))
# out_nlminb
# end_time <- Sys.time()
# end_time-start_time
# 
# 
# calibrated_nlminb=ParametersReconstruction(out_nlminb$par,n=4,common = FALSE)
# 
# correlation = cov2cor(calibrated_nlminb$S)
# correlation
# calibrated_nlminb$theta

calibrated_params = CalibrateMVMerton(x=cbind(eur[1:N],gbp[1:N],chf[1:N],jpy[1:N]), n=4, dt = dt )
#calibrated_params
cov2cor(calibrated_params$S)



###########################################################################
######################## Building complete correlaion #####################
###########################################################################
# from a 4 asset model


# ONLY WORKS IF N_ASSET IS *EVEN*
#N_assets = dim(my_returns)[2]/2
N_assets = 4


idx_matrix = matrix(seq(from = 1, to = N_assets, by = 1),N_assets%/% 2,2, byrow = TRUE)

final_cov = matrix(rep(0,N_assets*N_assets),N_assets,N_assets)

incremental = matrix(rep(0,N_assets*N_assets),N_assets,N_assets)
w_matrix = matrix(rep(0,N_assets*N_assets),N_assets,N_assets)

beg <- Sys.time()
for (i in 1:(N_assets%/% 2 -1)){
  for (j in (i+1):(N_assets%/% 2)){
    idx =c(idx_matrix[i,],idx_matrix[j,])
    
    w_matrix[idx,idx]=
      w_matrix[idx,idx]+1
    SS = CalibrateMVMerton(x=cbind(my_returns[,2*idx[1]],my_returns[,2*idx[2]],my_returns[,2*idx[3]],my_returns[,2*idx[4]]),
                           n=4, dt = dt )$S
    final_cov[idx,idx]= final_cov[idx,idx] + SS
    print(final_cov)
  }
}

end<- Sys.time()

w_matrix
final_cov/w_matrix

cov2cor(final_cov/w_matrix)*100

end-beg
