

library(mvtnorm)
library(DEoptim)
library(statmod)
library(NMOF)



library(xlsx)
library(openxlsx)


library(binaryLogic)

source('MultivariateMertonModel.R')
source('CalibrationMVMerton.R')


# CREATION OF THE DATASET OF LOG-RETURNS FROM EXCEL DATASET

# my_data<-read.xlsx("XBT_Correlations_nasdaq.xlsm",sheet = "DATA" , colNames =TRUE)
# 
# leng = dim(my_data)[1]
# N_assets = dim(my_data)[2] /2
# 
# var_names = c('btc_date','btc', #btc
#               'bric_date','bric','sp500_date','sp500','eurostoxx_date','eurostoxx','nasdaq_date','nasdaq', #stock
#               'bond_europe_date','bond_europe','bond_us_date','bond_us','bond_eur_date','bond_eur', #bond
#               'eur_date','eur','gbp_date','gbp','chf_date','chf','jpy_date','jpy', #FX
#               'gold_date','gold','wti_date','wti','grain_date','grain','metal_date','metal', # commodities
#               'vix_date','vix')
# 
# my_data = my_data[,var_names]
# 
# my_returns_values = log(my_data[1:(leng-1),2*(1:N_assets)]/my_data[2:leng,2*(1:N_assets)] )
# my_returns =my_data[1:(leng-1),]
# 
# my_returns[,2*(1:N_assets)]=my_returns_values
# 
# 
# save(my_returns, file = "returns.Rda")
# save(my_data, file = "data.Rda")

load("returns.Rda")
load("data.Rda")


N=500

# Plot of first few data
x11()
par(mfrow = c(2,3))
plot(my_data$btc_date[1:N],my_data$btc[1:N],type='l')
plot(my_data$nasdaq_date[1:N],my_data$nasdaq[1:N],type = 'l',col='green')
plot(my_data$eurostoxx_date[1:N],my_data$eurostoxx[1:N],type = 'l',col='blue')
plot(my_returns$btc_date[1:N],my_returns$btc[1:N], type='l')
plot(my_returns$nasdaq_date[1:N],my_returns$nasdaq[1:N], type='l',col='green')
plot(my_returns$eurostoxx_date[1:N],my_returns$eurostoxx[1:N], type='l',col='blue')

graphics.off()
  





################ TEST CALIBRATION on 4 assets ############################
attach(my_returns)


N = dim(my_returns)[1] #total number of observations

dt = 1/255




calibrated_params = CalibrateMVMerton(x=cbind(btc[1:N],bric[1:N],sp500[1:N],eurostoxx[1:N]), n=4, dt = dt )
calibrated_params

print(calibrated_params$mu)
colMeans(my_returns[,2*(1:4)]*255)

cov2cor(calibrated_params$S)


x11()
par(mfrow=c(2,2))
plot(my_data$btc_date, my_data$btc, type = 'l', main="BITCOIN", ylab = 'Price')
plot(my_data$bric_date, my_data$bric, type='l', main= "BRIC", ylab = 'Price')
plot(my_data$bond_us_date, my_data$bond_us, type='l', main= "BOND_US", ylab = 'Price')
plot(my_data$bond_eur_date, my_data$bond_eur, type='l', main= "BOND_EUR", ylab = 'Price')

x11()
pairs(my_returns[,(1:4)*2])




# x11()
# par(mfrow=c(2,2))
# for (idx in 1:4) {
#   xx=seq(-2,2,length.out = 2000)
#   yy= MultivariateMertonPdf_1asset_nocommon(x = xx,dt = 1/255,mu = calibrated_params$mu[idx], S = calibrated_params$sigma[idx],
#                                         theta = calibrated_params$theta[idx], delta = calibrated_params$delta[idx], lambda = calibrated_params$lambda[idx])
#   
#   print(sum(yy*(xx[2]-xx[1])))
#   print(sum(xx*yy*(xx[2]-xx[1]))*255)
#   
#   hist(my_returns[,idx*2], breaks = 60,freq = FALSE)
#   lines(xx, dnorm(xx, mean = mean(my_returns[,idx*2]),sd = sd(my_returns[,idx*2])))
#   lines(xx,yy, col='blue')
# }



###########################################################################
######################## Building complete correlation #####################
###########################################################################
# from a 4 asset model


N_assets = dim(my_returns)[2]/2

initial_time = Sys.time()

# First work with an even number of assets
N_assets_even = ifelse(N_assets %%2 ==1, N_assets -1, N_assets)

N_assets_even

N=dim(my_returns)[1] 

idx_matrix = matrix(seq(from = 1, to = N_assets_even, by = 1),N_assets_even%/% 2,2, byrow = TRUE)

final_cov = matrix(rep(0,N_assets_even*N_assets_even),N_assets_even,N_assets_even)


w_matrix = matrix(rep(0,N_assets_even*N_assets_even),N_assets_even,N_assets_even)

beg <- Sys.time()



N_computations = (N_assets_even %/% 2)*(N_assets_even %/% 2 -1) /2

names = colnames(my_returns[,2*(1:N_assets_even)])
mus = matrix(rep(0,N_computations*N_assets_even),N_computations,N_assets_even)
colnames(mus) = names
thetas = matrix(rep(0,N_computations*N_assets_even),N_computations,N_assets_even)
colnames(thetas) = names
deltas = matrix(rep(0,N_computations*N_assets_even),N_computations,N_assets_even)
colnames(deltas) = names
lambdas = matrix(rep(0,N_computations*N_assets_even),N_computations,N_assets_even)

colnames(lambdas) = names

ctr = 1
for (i in 1:(N_assets_even%/% 2 -1)){
  for (j in (i+1):(N_assets_even%/% 2)){
    idx =c(idx_matrix[i,],idx_matrix[j,])
    
    w_matrix[idx,idx]= w_matrix[idx,idx]+1
    calibrated = CalibrateMVMerton(x=cbind(my_returns[,2*idx[1]][1:N],my_returns[,2*idx[2]][1:N],my_returns[,2*idx[3]][1:N],my_returns[,2*idx[4]][1:N]),
                                   n=4, dt = dt )
    
    mus[ctr,idx] = calibrated$mu
    thetas[ctr,idx] = calibrated$theta
    deltas[ctr,idx] = calibrated$delta
    lambdas[ctr,idx] = calibrated$lambda
    ctr = ctr+1
    
    mus[ctr,idx] = calibrated$mu
    thetas[ctr,idx] = calibrated$theta
    deltas[ctr,idx] = calibrated$delta
    lambdas[ctr,idx] = calibrated$lambda
    ctr = ctr+1
    
    SS = calibrated$S
    final_cov[idx,idx]= final_cov[idx,idx] + SS
    # print(final_cov)
    out = capture.output(calibrated)
    cat(paste("Assets:", idx,sep = " "), out, 
        file = paste0("computation_of_full_corr_matrix_nasdaq",format(Sys.time(), "%Y-%m-%d"),".txt"),
        sep="\n", append=TRUE)
    }
}

# # Visual check of results
# w_matrix
# final_cov/w_matrix
covariance = final_cov/w_matrix


# ********* If total assets are ODD, need to calibrated the remaining one ***********
# N_assets = dim(my_returns)[2]/2

if(N_assets %% 2 == 1){

  complete_var = matrix(rep(0, N_assets*N_assets), nrow = N_assets)
  complete_var[1:(N_assets-1), 1:(N_assets-1)] = covariance * w_matrix
  
  complete_w_matrix = matrix(rep(0, N_assets*N_assets), nrow = N_assets)
  complete_w_matrix[1:(N_assets-1), 1:(N_assets-1)] = w_matrix
  
  idx_last_column = matrix(1:(3*(N_assets%/%3)), ncol = 3, byrow = TRUE)
  idx_last_column = rbind(idx_last_column, N_assets - c(3,2,1))
  
  last_mu=0
  last_theta = 0
  last_delta = 0
  last_lambda= 0 
  for (i in 1:dim(idx_last_column)[1]){
    
    idx =c(idx_last_column[i,],N_assets)
    print(idx)
    complete_w_matrix[idx,idx] = complete_w_matrix[idx,idx]+1
    
    calibrated = CalibrateMVMerton(x=cbind(my_returns[,2*idx[1]][1:N],my_returns[,2*idx[2]][1:N],my_returns[,2*idx[3]][1:N],my_returns[,2*idx[4]][1:N]),
                                   n=4, dt = dt )

    last_mu = last_mu+  calibrated$mu[4]
    last_theta = last_theta+ calibrated$theta[4]
    last_delta = last_delta+ calibrated$delta[4]
    last_lambda = last_lambda+ calibrated$lambda[4]
    
    SS = calibrated$S
    complete_var[idx,idx]= complete_var[idx,idx] + SS
    # print(final_cov)
    out = capture.output(calibrated)
    cat(paste("Assets:", idx,sep = " "), out, 
        file = paste0("computation_of_full_corr_matrix_nasdaq",format(Sys.time(), "%Y-%m-%d"),".txt"),
        sep="\n", append=TRUE)
  }
  
  
  last_mu = last_mu/dim(idx_last_column)[1]
  last_theta = last_theta/dim(idx_last_column)[1]
  last_delta = last_delta/dim(idx_last_column)[1]
  last_lambda = last_lambda/dim(idx_last_column)[1]
  
  complete_w_matrix
  covariance = complete_var/complete_w_matrix
  
  correlation = cov2cor(covariance)
  
  out = capture.output(cov2cor(final_cov/w_matrix))
  cat("FINAL RESULT FOR CORRELATION", out, 
      file = paste0("computation_of_full_corr_matrix_nasdaq",format(Sys.time(), "%Y-%m-%d"),".txt"),
      sep="\n", append=TRUE)
}



end<- Sys.time()
end-beg


results = list(covariance = covariance, correlation = correlation, 
               full_mu = mus, full_theta = thetas, full_delta= deltas, full_lambda = lambdas)

final_time = Sys.time()

final_time - initial_time

save(results, file= "results.Rda")


####### Save results to a .txt file #################
# 
# write.table(results$parameters, file = "results_txt_merton.txt")
# write.table(results$covariance, file = "results_txt_merton.txt", append = TRUE)
# write.table(results$correlation, file = "results_txt_merton.txt", append = TRUE)
# write.table(results$mus, file = "results_txt_merton.txt", append = TRUE)
# write.table(results$thetas, file = "results_txt_merton.txt", append = TRUE)
# write.table(results$deltas, file = "results_txt_merton.txt", append = TRUE)



