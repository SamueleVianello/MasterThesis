###untitled####
## # # # # # # # 

source("MertonParsimoniousCalibration.R")
source("MultivariateMertonModel.R")
source("CalibrationMVMerton.R")
load("returns.Rda")
attach(my_returns)


N = dim(my_returns)[1] #total number of observations
N_assets= 16
dt = 1/255


full_results = matrix(0, nrow = N_assets, ncol = 5)
rownames(full_results) = colnames(my_returns)[2*(1:N_assets)]
colnames(full_results) = c("mu", "sigma","theta","delta","lambda")
pdf_mean = matrix(0, nrow = N_assets, ncol = 1)
rownames(pdf_mean) = colnames(my_returns)[2*(1:N_assets)]
  
for (i in 1:N_assets){
  asset = my_returns[,2*i]
  asset_name = colnames(my_returns)[2*i]
  
  calibrated_params = CalibrateMVMerton(x=matrix(asset[1:N], ncol=1), n=1, dt = dt, custom_jump_bounds = T)
  
  params = matrix(c(calibrated_params$mu-calibrated_params$S/2, sqrt(calibrated_params$S), calibrated_params$theta, calibrated_params$delta, calibrated_params$lambda), nrow = 1)
  params = rbind( params, c(mean(asset)*255, sd(asset)*sqrt(255), 0,0,0))
  colnames(params) = c("mu-0.5*sigma^2", "sigma","theta","delta","lambda")
  params

  full_results[i,]=c(calibrated_params$mu, sqrt(calibrated_params$S), calibrated_params$theta, calibrated_params$delta, calibrated_params$lambda)

  xx = seq(-0.6,0.6, length.out = 4000)
  yy_merton = MultivariateMertonPdf_1asset_nocommon(x = xx,dt = dt, mu = calibrated_params$mu,S = calibrated_params$S,
                                                    theta = calibrated_params$theta, delta = calibrated_params$delta, lambda = calibrated_params$lambda)
  yy_gauss = dnorm(xx, mean = mean(asset), sd = sd(asset))

  hist(asset, freq = F, breaks = 40 , main = asset_name )
  grid()
  lines(xx,yy_merton,col='blue', type = 'l')
  lines(xx,yy_gauss)

  res = matrix(c(sum(xx*yy_merton*(xx[2]-xx[1])), mean(asset)), ncol=1)

  abline(v=sum(xx*yy_merton*(xx[2]-xx[1])), col = 'blue', lwd=3)
  abline(v=mean(asset))

  pdf_mean[i,1]=sum(xx*yy_merton*(xx[2]-xx[1]))/dt
  cbind(res, res*255)
  params
}
  
full_results

# calibrating model correlation matrix

n=16

corr_matrix = cor(my_returns[,2*(1:n)])
beg = Sys.time()
c_mat = calibrate_full_correlation_merton(sample_corr_matrix = corr_matrix, Nasset = n,
                                  mu =full_results[1:n,1] , vol = full_results[1:n,2],
                                  mu_j = full_results[1:n,3], sigma_j = full_results[1:n,4], lambda = full_results[1:n,5],
                                  dt = 1/255,final_t = 2, Nsim = 1000)
Sys.time() -beg
c_mat
corr_matrix[1:n,1:n]



  
