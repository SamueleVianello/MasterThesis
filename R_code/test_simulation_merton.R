load("MertonParsimoniousCalibration.R")

calibrated_params = CalibrateMVMerton(x=cbind(nasdaq,sp500), n=2, dt = dt, custom_jump_bounds = T)
jj=2

# one-asset simulation
sim_x = simulate_mv_merton(Nasset = 1,x0 = 0, mu =calibrated_params$mu[jj] ,vol = sqrt(calibrated_params$S[jj,jj]),
                   mu_j = calibrated_params$theta[jj], sigma_j =  calibrated_params$delta[jj], lambda = calibrated_params$lambda[jj],
                  Nsim = 1000)
# 2-asset simulation
sim_x2 = simulate_mv_merton(x0 = rep(0,length(calibrated_params$mu)), mu =calibrated_params$mu ,vol = sqrt(diag(calibrated_params$S)), CorrMatrix = calibrated_params$corr,
                           mu_j = calibrated_params$theta, sigma_j =  calibrated_params$delta, lambda = calibrated_params$lambda,
                           Nsim = 1000, final_t = 1)


matplot(t(sim_x$x_1[1:100,]), type = 'l')

# check that simulated returns follow the correct distribution
jj=2
hist(sim_x2[[jj]][,2], freq = FALSE, breaks = 40)
grid()
xx = seq(-0.6,0.6, length.out = 1000)
yy_merton = MultivariateMertonPdf_1asset_nocommon(x = xx,dt = dt, mu =calibrated_params$mu[jj] ,S = calibrated_params$S[jj,jj],
                                                    theta = calibrated_params$theta[jj], delta =calibrated_params$delta[jj], 
                                                  lambda = calibrated_params$lambda[jj])
lines(xx,yy_merton,col='blue')


cor_est = correlation_MC_estimation(sim_x2)
cor_est
(cor(cbind(sp500, btc, bric, metal))-cor_est)/cor(cbind(sp500, btc, bric, metal))



# estimate correlation from model parameters
expected_model_correlation(model_corr=0.5284209,mu =calibrated_params$mu ,vol = sqrt(diag(calibrated_params$S)),
                            mu_j = calibrated_params$theta, sigma_j = calibrated_params$delta, lambda = calibrated_params$lambda,
                                      x0=rep(0,2), dt=1/255, final_t=1, Nsim=1000)
# OK


# calibrate model correlation from sample corr [2 assets]
calibrate_correlation_merton(sample_corr = 0.528, Nasset = 2, mu =calibrated_params$mu, vol = sqrt(diag(calibrated_params$S)),
                             mu_j = calibrated_params$theta, sigma_j = calibrated_params$delta, lambda = calibrated_params$lambda,
                             dt = 1/255,final_t = 1,Nsim = 1000)
# OK

expected_model_correlation(model_corr = -1,x0=rep(0,2),#Nasset = n-1,
mu =full_results[2:n,1] , vol = full_results[2:n,2],
mu_j = full_results[2:n,3], sigma_j = full_results[2:n,4], lambda = full_results[2:n,5],
dt = 1/255, final_t = 1,Nsim = 1000)


# calibrate using full function
calibrate_full_correlation_merton(sample_corr_matrix = matrix(c(1,0.3934197,0.3934197,1), ncol=2),Nasset = 2,
                                  mu =calibrated_params$mu, vol = sqrt(diag(calibrated_params$S)),
                                  mu_j = calibrated_params$theta, sigma_j = calibrated_params$delta, lambda = calibrated_params$lambda,
                                  dt = 1/255,final_t = 1,Nsim = 1000)



#test full correlation calibration
n=16

corr_matrix = cor(my_returns[,2*(1:n)])
beg = Sys.time()
calibrate_full_correlation_merton(sample_corr_matrix = corr_matrix, Nasset = n,
                                  mu =full_results[1:n,1] , vol = full_results[1:n,2],
                                  mu_j = full_results[1:n,3], sigma_j = full_results[1:n,4], lambda = full_results[1:n,5],
                                  dt = 1/255,final_t = 1,Nsim = 1000)
corr_matrix
Sys.time() -beg
