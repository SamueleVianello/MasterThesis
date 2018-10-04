###########################################
##### Calibration of Bates and Heston #####
###########################################

# need cumulative returns X_t =  log (S_t/S_0)

CalibrateHeston=function(x, x_0,sigma_0, dt, trace = 10){
  
  library(DEoptim)
  source("BatesModel.R")
  
  bounds = BoundsCreator()
  
  # First optimization using deoptim
  print("Starting calibration using DEoptim...")
  control_list_deoptim = list(itermax = 50, NP = 200, strategy = 1,trace=trace)
  
  start_time_deoptim <- Sys.time()
  outDE <- DEoptim(negloglikHeston, 
                   lower = bounds$lower, upper = bounds$upper, control = control_list_deoptim, 
                    x=x, x_0 = x_0, sigma_0=sigma_0, dt = dt)
  end_time_deoptim <- Sys.time()
  
  # Second and final optimization using nlminb
  print("Starting calibration using nlminb...")
  
  initial=outDE$optim$bestmem
  start_time_nlminb <- Sys.time()
  out_nlminb = nlminb(initial,objective = negloglikHeston,lower = bounds$lower, upper = bounds$upper,
                      x=x, x_0 = x_0, sigma_0=sigma_0,dt = dt,
                      control=list(eval.max = 10000,iter.max = 1000, trace = trace))
  print(out_nlminb)
  end_time_nlminb <- Sys.time()
  
  print(paste("DEoptim:",end_time_deoptim - start_time_deoptim))
  print(paste("nlminb:", end_time_nlminb - start_time_nlminb))
  print(paste("TOTAL:", end_time_nlminb - start_time_deoptim))
  
  res = ParametersReconstruction(out_nlminb$par)
  res[["message"]] = out_nlminb$message
  res[["objective_function"]] = out_nlminb$objective
  res[["total_time"]] = end_time_nlminb - start_time_deoptim
  return(res)
}


BoundsCreator= function(n=1){
  # Creates lower and upper boundaries for the DEoptim optimization on the likelihood
  # for a n-multivariate merton process and n_common common jumps
  min_mu = -5
  max_mu = 5
  # min_lambda = 0.1
  # max_lambda = 100
  min_sigma = 1e-5
  max_sigma = 10
  min_corr = -1
  max_corr = 1

  low = c(min_mu,min_mu,min_sigma,min_sigma,min_corr)
  up = c(max_mu,max_mu,max_sigma,max_sigma,max_corr)
  return(list(lower = low, upper = up))
}

ParametersReconstruction = function(params){
  return(list(r=params[1],k=params[2],eta=params[3],theta=params[4],rho=params[5]))
}






