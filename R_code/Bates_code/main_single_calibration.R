library(PerformanceAnalytics)
source("BatesModel.R")
source("AuxiliaryFunctions.R")
source("CalibrationBates.R")
load("returns.Rda")


attach(my_returns)

N = dim(my_returns)[1]
dn = N  # in case you don't want to use the whole time series, change it



asset = chf[1:dn]
asset_date = chf_date[1:dn]
asset_name = 'chf'

t_very_beginning = Sys.time()

# Uncomment if calibrating on multiple assets
asset_idx = 5:14
result_H = matrix(ncol=10, nrow = length(asset_idx))
result_B = matrix(ncol=13, nrow = length(asset_idx))
# begin of for loop to cycle through all assets
for(j in asset_idx){ #**********************
asset = my_returns[1:dn,2*j]
asset_date = my_returns[1:dn, 2*j-1]
asset_name = colnames(my_returns)[2*j]




# General plot
plot_returns(asset_date[1:dn],asset[1:dn], asset_name)






# Setting initial values [required only if global optimization is not included]
initial_mu = mean(asset[1:dn])*255
initial_var = var(asset[1:dn])*255

initial_mu_j = unname(quantile(asset[1:dn], probs = 0.001))
initial_lambda = sum(asset[1:dn] <= initial_mu_j)*255/dn

initial = c(r=initial_mu, k= 1.5, eta = initial_var*0.1, theta = 0.1, rho = -0.3, mu_j =drop(initial_mu_j), sigma_j=0.02, lambda=initial_lambda)
# Making sure feller condition is verified, otherwise modifying k and theta as suggested by the warning
initial = check_feller(initial)



# Creating range for rho and  mu_jump based on skewness of returns

alpha = 0.999 # quantile for max jump

sample_skew = skewness(asset[1:dn])*sqrt(dn*(dn-1))/(dn-2)
SES = sqrt(6*dn*(dn-1)/ ((dn-2)*(dn+1)*(dn+3)))
skew_significance = sample_skew/SES

if(skew_significance <= -2) {
  min_muj = 2*quantile(asset[1:dn], probs = 1-alpha)
  max_muj = quantile(asset[1:dn], probs = 1-alpha)
  min_rho = -1
  max_rho = -0.0001
} else if(skew_significance >-2){
  min_muj = quantile(asset[1:dn], probs = alpha)
  max_muj = 2*quantile(asset[1:dn], probs = alpha)
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

# Saving output to file with asset name
write.table(x = res_H,file = paste0("params_", asset_name, ".txt"), sep="\t")
write.table(x = res_B,file = paste0("params_", asset_name, ".txt"), sep="\t",append = TRUE)

result_B[j,] = res_B[1,]     # Add this only if cycling on multiple assets
} #*********** end for to cycle through all assets

View(res_H)
View(res_B)


t_very_end = Sys.time()

# total computation time for everything
t_very_end - t_very_beginning










# Plot of results as a visual comparison ------------------------------------------------------
par_heston = list_from_array(colMeans(res_H)[1:5], model='heston')
par_bates = list_from_array(colMeans(res_B)[1:8], model='bates')


hist(asset[1:dn], freq = FALSE, breaks = 200, main=paste("Histogram of", asset_name))

x_max = 2*max(abs(asset))
xx=seq(-x_max,x_max, length.out = 4000)

yy_H = pdfHeston(x=xx, dt=(1/255), r=par_heston$r, k=par_heston$k, eta=par_heston$eta, theta=par_heston$theta, rho=par_heston$rho,N=2**12)
lines(xx,yy_H,col='red',type='l')
sum(yy_H*(xx[2]-xx[1]))

yy = pdfBates(x=xx, dt=(1/255), r=par_bates$r, k=par_bates$k, eta=par_bates$eta, theta=par_bates$theta, rho=par_bates$rho, 
                     mu_j = par_bates$mu_j, sigma_j = par_bates$sigma_j, lambda = par_bates$lambda,N=2**12 )
lines(xx,yy,col='blue',type='l')
sum(yy*(xx[2]-xx[1]))

# Add gaussian fit
lines(xx, dnorm(xx,mean = mean(asset[1:dn]),sd = sd(asset[1:dn])))
abline(v=mean(asset), lty=4)
legend("topleft", legend = c('gaussian', 'heston', 'bates'), col=c('black', 'red', 'blue'), lty=c(1,1,1))



# compute relevant moments

sample_stat = c(sample_mean = mean(asset[1:dn]),sample_variance = var(asset[1:dn]),
                sample_skewness=skewness(asset[1:dn]),sample_kurtosis=kurtosis(asset[1:dn]) )

heston_stat = c(model_mean_heston = mean_from_pdf(xx,yy_H),model_variance_heston=variance_from_pdf(xx,yy_H),
                model_skewness_heston = skewness_from_pdf(xx,yy_H),model_kurtosis_heston = kurtosis_from_pdf(xx,yy_H))

bates_stat= c(model_mean_bates=mean_from_pdf(xx,yy),model_variance_bates = variance_from_pdf(xx,yy),
              model_skewness_bates = skewness_from_pdf(xx,yy),model_kurtosis_bates = kurtosis_from_pdf(xx,yy))

sample_stat
heston_stat
bates_stat








# Computing empirical moments for all assets ------------------------------------------------

N_assets = (dim(my_returns)[2]/2)
dn = N

moments_matrix = matrix(rep(0,4*N_assets), ncol = 4 )

for(i in 1: N_assets){
  asset = my_returns[1:dn,2*i]
  asset_name = colnames(my_returns)[2*i]
  sample_stat = c(sample_mean = mean(asset[1:dn]),sample_variance = var(asset[1:dn]),
                  sample_skewness=skewness(asset[1:dn]),sample_kurtosis=kurtosis(asset[1:dn]) )
  moments_matrix[i,]=sample_stat
}
rownames(moments_matrix)= colnames(my_returns[,2*(1:N_assets)])
View(moments_matrix)
  
  

