library(PerformanceAnalytics)
source("BatesModel.R")
source("AuxiliaryFunctions.R")
source("CalibrationBates.R")
source("SimulationBates.R")
load("returns.Rda")




attach(my_returns)

N = dim(my_returns)[1] 
dn = N  # in case you don't want to use the whole time series, change this

N_assets = 16
result_H = matrix(ncol=10, nrow = N_assets)
result_B = matrix(ncol=13, nrow = N_assets)



# begin of for loop to cycle through all assets
for(j in 13:16){ #**********************
  asset = my_returns[1:dn,2*j]
  asset_date = my_returns[1:dn, 2*j-1]
  asset_name = colnames(my_returns)[2*j]



  # General plot
  # plot_returns(asset_date[1:dn],asset[1:dn], asset_name)




  # Setting initial values [required only if global optimization is not included]
  initial_mu = mean(asset[1:dn])*255
  initial_var = var(asset[1:dn])*255

  initial_mu_j = unname(quantile(asset[1:dn], probs = 0.001))
  initial_lambda = sum(asset[1:dn] <= initial_mu_j)*255/dn

  initial = c(r=initial_mu, k= 1.5, eta = initial_var*0.1, theta = 0.1, rho = -0.3, mu_j =drop(initial_mu_j), sigma_j=0.02, lambda=initial_lambda)
  # Making sure feller condition is verified, otherwise modifying k and theta as suggested by the warning
  initial = check_feller(initial)



  # Creating range for rho and  mu_jump based on skewness of returns
  alpha_max = 0.995 # quantile for max jump
  alpha_min = 0.999 # quantile for min jump

  sample_skew = skewness(asset[1:dn])*sqrt(dn*(dn-1))/(dn-2)
  SES = sqrt(6*dn*(dn-1)/ ((dn-2)*(dn+1)*(dn+3)))
  skew_significance = sample_skew/SES

  if(skew_significance <= -2) {
    min_muj = 2*quantile(asset[1:dn], probs = 1-alpha_min)
    max_muj = quantile(asset[1:dn], probs = 1-alpha_max)
    min_rho = -1
    max_rho = -0.0001
  } else if(skew_significance >-2){
    min_muj = quantile(asset[1:dn], probs = alpha_max)
    max_muj = 2*quantile(asset[1:dn], probs = alpha_min)
    min_rho = 0.0001
    max_rho = +1
  } else{
    min_muj = 2*min(asset[1:dn])       # choose what is a significant bound for jumps
    max_muj = 2*max(asset[1:dn])       # choose what is a significant bound for jumps
    min_rho = -1
    max_rho = +1
  }


  # Calibrating multiple times to check if results are stable ----------------------------------------

  Nrep = 1   # modify to change number of repetitions

  res_H = matrix(nrow = Nrep, ncol=10)  # 5 parameters + 1 NegLogLik + 4 moments = 10
  colnames(res_H) = c('r','k','eta','theta','rho', 'NegLogLik', 'mean', 'var', 'skew', 'kurt')
  res_B = matrix(nrow = Nrep, ncol=13)  # 8 parameters + 1 NegLogLik + 4 moments = 13
  colnames(res_B) = c('r','k','eta','theta','rho', 'mu_j','sigma_j','lambda','NegLogLik', 'mean', 'var', 'skew', 'kurt')

  t_start = Sys.time()
  # ____________________________Heston calibration__________________________________
  for(i in 1:Nrep){
    t1 = Sys.time()
    par_heston = CalibrateModel(x=asset[1:dn], dt=1/255, trace=10,
                                initial = c(initial[1], initial[2:5]),
                                global = T,local = T, feller = T , model="heston",
                                fixing = c(1),
                                set_upper = 5, upper_bnd = max_rho,
                                set_lower = 5, lower_bnd = min_rho)
    t2=Sys.time()
    t2-t1

    xlim = max(abs(asset[1:dn]))*2
    xx=seq(-xlim, xlim, length.out = 4000)

    yy_H = pdfHeston(x=xx, dt=(1/255), r=par_heston$r, k=par_heston$k, eta=par_heston$eta, theta=par_heston$theta, rho=par_heston$rho,N=2**12)
    sum(yy_H*(xx[2]-xx[1]))

    heston_stat = c(model_mean_heston = mean_from_pdf(xx,yy_H),model_variance_heston=variance_from_pdf(xx,yy_H),
                    model_skewness_heston = skewness_from_pdf(xx,yy_H),model_kurtosis_heston = kurtosis_from_pdf(xx,yy_H))

    res_H[i,]=array(as.numeric(unlist(c(par_heston, heston_stat))), dim = c(1,10))

  }

  result_H[j,] = res_H[1,]     # Add this only if cycling on multiple assets

  t_mid = Sys.time()
  print("***********************************")
  print("End of calibration of Heston Model.")
  print("***********************************")

  # _________________________Bates Calibration______________________________________
  for(i in 1:Nrep){
    t1 = Sys.time()
    min
    par_bates = CalibrateModel(x=asset[1:dn], dt=1/255, trace=10,
                               initial = (c(initial[1], initial[2:8])),
                               global =T, local = T,  feller = T ,model='bates',
                               fixing = 1,
                               set_upper=c(5,6), upper_bnd=c(max_rho, max_muj),
                               set_lower = c(5,6), lower_bnd = c(min_rho, min_muj))
    t2=Sys.time()
    t2-t1

    xlim = max(abs(asset[1:dn]))*2
    xx=seq(-xlim, xlim, length.out = 4000)
    yy = pdfBates(x=xx, dt=(1/255), r=par_bates$r, k=par_bates$k, eta=par_bates$eta, theta=par_bates$theta, rho=par_bates$rho,
                  mu_j = par_bates$mu_j, sigma_j = par_bates$sigma_j, lambda = par_bates$lambda,N=2**12 )
    sum(yy*(xx[2]-xx[1]))
    bates_stat= c(model_mean_bates=mean_from_pdf(xx,yy),model_variance_bates = variance_from_pdf(xx,yy),
                  model_skewness_bates = skewness_from_pdf(xx,yy),model_kurtosis_bates = kurtosis_from_pdf(xx,yy))
    res_B[i,]=array(as.numeric(unlist(c(par_bates, bates_stat))), dim = c(1,13))
  }

  t_end = Sys.time()

  print("************** total time **************")
  print('Heston:')
  print(t_mid - t_start)
  print('Bates:')
  print(t_end - t_mid)
  print('Total:')
  print(t_end - t_start)


  result_B[j,] = res_B[1,]     # Add this only if cycling on multiple assets

  x11()
  hist(asset, freq = FALSE, breaks = 60, main = asset_name)
  lines(xx,yy_H, col='red')
  lines(xx,yy, type = 'l',col='blue' ) 
  abline(v=sum(xx*yy_H*(xx[2]-xx[1])), col='red')
  abline(v=sum(xx*yy*(xx[2]-xx[1])), col='blue')
} #*********** end for to cycle through all assets



result_H
result_B


# save(file = 'single_calibration_Heston_Bates.Rda', list = c('result_B', 'result_H'))





# ____________________ CALIBRATING MODEL CORRELATION __________________________________

N_assets=16
# computing sample correlation from data
sample_corr  = cor(my_returns[,2*(1:N_assets)])
sample_corr

# getting an estimate of initial variance
idx = N:(N-10) 
var0 = diag(var(my_returns[idx,2*(1:N_assets)]))*255

# calibrating model correlation matrix
model_corr_H = calibrate_full_correlation(sample_corr_matrix = sample_corr, Nasset = N_assets,
                                        mu = result_H[1:N_assets,1], k=result_H[1:N_assets,2], eta=result_H[1:N_assets,3], theta=result_H[1:N_assets,4], rho=result_H[1:N_assets,5],
                                        V0 = var0, dt=1/255, final_t = 1, Nsim = 400)


# simulating path 
mvsim_H = simulate_mv_heston(Nasset = N_assets,
                           mu = result_H[1:N_assets,1], k=result_H[1:N_assets,2], eta=result_H[1:N_assets,3], theta=result_H[1:N_assets,4], rho=result_H[1:N_assets,5],
                           CorrMatrix = model_corr_H,
                           V0 = var0, dt=1/255, final_t = 1, Nsim = 400)

# check on correlation from simulated path and original sample correlation 
correlation_MC_estimation(mvsim_H)
sample_corr



# Same thing for Bates
# calibrating model correlation matrix
model_corr_B= calibrate_full_correlation(sample_corr_matrix = sample_corr, Nasset = N_assets,
                                          mu = result_B[1:N_assets,1], k=result_B[1:N_assets,2], eta=result_B[1:N_assets,3], theta=result_B[1:N_assets,4], rho=result_B[1:N_assets,5],
                                          mu_j = result_B[1:N_assets,6], sigma_j = result_B[1:N_assets,7], lambda = result_B[1:N_assets,8],
                                          V0 = var0, dt=1/255, final_t = 2, Nsim = 500)
model_corr_B


# simulating path 
mvsim_B = simulate_mv_bates(Nasset = N_assets,
                             mu = result_B[1:N_assets,1], k=result_B[1:N_assets,2], eta=result_B[1:N_assets,3], theta=result_B[1:N_assets,4], rho=result_B[1:N_assets,5],
                             mu_j = result_B[1:N_assets,6], sigma_j = result_B[1:N_assets,7], lambda = result_B[1:N_assets,8],
                             CorrMatrix = model_corr_B,
                             V0 = var0, dt=1/255, final_t = 1, Nsim = 500)

# check on correlation from simulated path and original sample correlation 
correlation_MC_estimation(mvsim_B)
sample_corr




save(file = 'calibration_Heston_Bates.Rda', list = c('result_H','model_corr_H', 'result_B', 'model_corr_B' ))




# BONUS ***
# mean reversion level (eta) estimator:
etas= colMeans(my_returns[,2*(1:14)]^2)*255



