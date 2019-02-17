###############################################################
############ Calibration of a Multivariate Merton #############
###############################################################


CalibrateMVMerton=function(x, n, dt, trace = 10, custom_jump_bounds =T){
  
  library(DEoptim)
  source("MultivariateMertonModel.R")
  
  # Choose objective function given the number of assets
  if(n ==1) obj = negloglik_1asset_nocommon
  else if(n ==2) obj = negloglik_2assets_nocommon
  else if(n ==3) obj = negloglik_3assets_nocommon
  else if(n ==4) obj = negloglik_4assets_nocommon
  else obj = negloglik
  
  if(custom_jump_bounds){
    min_jump = rep(0,n)
    max_jump = rep(0,n)
    
    alpha_max = 0.995 # quantile for max jump
    alpha_min = 0.999 # quantile for min jump
    for (i in 1:n) {
      min_jump[i] = 2*quantile(x= x[,i], probs = 1-alpha_min)
      max_jump[i] = quantile(x= x[,i], probs = 1-alpha_max)
    }
  }
  else{
    min_jump = rep(-0.1,n) # default jumps at -10%
    max_jump = rep(-1,n)
  }
  
  bounds_nocommon = BoundsCreator(n=n, custom_jump_mean = custom_jump_bounds, 
                                  max_jump_mean = max_jump, min_jump_mean = min_jump)
  
  print(rbind(bounds_nocommon$lower, bounds_nocommon$upper))
  
  # First optimization using deoptim
  print("Starting calibration using DEoptim...")
  control_list_deoptim = list(itermax = 50, NP = 10*length(bounds_nocommon$lower), strategy = 6,trace=trace)
  
  start_time_deoptim <- Sys.time()
  outDE <- DEoptim(obj, 
                   lower = bounds_nocommon$lower, upper = bounds_nocommon$upper, control = control_list_deoptim, 
                   dt = dt, x = x, n = n)
  end_time_deoptim <- Sys.time()
  
  # Second and final optimization using nlminb
  print("Starting calibration using nlminb...")
  
  initial=outDE$optim$bestmem
  start_time_nlminb <- Sys.time()
  out_nlminb = nlminb(initial,objective = obj,lower = bounds_nocommon$lower,
                      upper = bounds_nocommon$upper,dt=dt, x=x, n=n,
                      control=list(eval.max = 10000,iter.max = 1000, trace = trace))
  print(out_nlminb)
  end_time_nlminb <- Sys.time()
  
  print(paste("DEoptim:",end_time_deoptim - start_time_deoptim))
  print(paste("nlminb:", end_time_nlminb - start_time_nlminb))
  print(paste("TOTAL:", end_time_nlminb - start_time_deoptim))
  
  res = ParametersReconstruction(out_nlminb$par,n=n,common = FALSE)
  res[["message"]] = out_nlminb$message
  res[["objective_function"]] = out_nlminb$objective
  res[["total_time"]] = end_time_nlminb - start_time_deoptim
  return(res)
}









BoundsCreator= function(n, 
                        custom_jump_mean =FALSE,
                        max_jump_mean = rep(NA,n),
                        min_jump_mean = rep(NA,n)
){
  # Creates lower and upper boundaries for the DEoptim optimization on the likelihood
  # for a n-multivariate merton process and n_common common jumps
  min_mu = -5
  max_mu = 5
  min_sigma = 1e-5
  max_sigma = 5
  min_lambda = 0.0001
  max_lambda = 100
  min_theta = -1
  max_theta = -0.1
  min_delta = 0.0001
  max_delta = 0.4
  min_corr = -1
  max_corr = 1
  
  
  if (n==1){
    low = c(min_mu, min_sigma, min_lambda, min_theta, min_delta)
    up = c(max_mu, max_sigma, max_lambda, max_theta, max_delta)
    
    if (custom_jump_mean){
      low[4] = min_jump_mean
      up[4] = max_jump_mean
    }
  }
  else {
    # initialising resulting vector low and up
    leng = 4*n + (n+1)*n*0.5
    low = rep(0,leng)
    up = rep(0,leng)
    
    
    idx =1
    # mean of continuos part
    low[idx:(idx+n-1)] = rep(min_mu,n)
    up[idx:(idx+n-1)] = rep(max_mu,n)
    idx = idx+n 
    
    # diffusion coefficients
    low[idx:(idx+n-1)] = rep(min_sigma,n)
    up[idx:(idx+n-1)] = rep(max_sigma,n)
    idx = idx+n 
    
    # correlations 
    N_var = n*(n-1)/2
    low[idx:(idx + N_var -1)] = rep(min_corr,N_var)
    up[idx:(idx + N_var -1)] = rep(max_corr,N_var)
    idx = idx + N_var
    
    # means of idiosyncratic jump term
    if (!custom_jump_mean){
      low[idx:(idx+n-1)] = rep(min_theta,n)
      up[idx:(idx+n-1)] = rep(max_theta,n)
    }
    else {
      for (i in 0:(n-1)) {
        low[idx+i] = ifelse(is.na(min_jump_mean[i+1]), min_theta, min_jump_mean[i+1])
        up[idx+i] = ifelse(is.na(max_jump_mean[i+1]), min_theta, max_jump_mean[i+1])
      }
    }
    idx = idx+n
    
    # standar deviation of idyosincratic jump term
    low[idx:(idx+n-1)] = rep(min_delta,n)
    up[idx:(idx+n-1)] = rep(max_delta,n)
    idx = idx+n
    
    # lambda of idyosincratic poissons
    low[idx:(idx+n-1)] = rep(min_lambda,n)
    up[idx:(idx+n-1)] = rep(max_lambda,n)
    idx = idx+n
  }
  
  if( n!=1 && ((idx-1)!=length(low))  || (length(low)!=length(up)) )
    stop("Error in parameter reconstruction: number of parameters is wrong.")
  
  return(list(lower = low, upper = up))
}



ParametersReconstruction = function(params, n, common = FALSE){
  
  if (n==1){
    mu = params[1]
    sigma = matrix(params[2])
    lambda = params[3]
    theta = params[4]
    delta = params[5]
    S=sigma^2
    corr=NA
  }
  else{
    # reconstruction of parameters:
    idx =1
    mu = params[idx:(idx+n-1)]
    idx = idx+n 
    
    
    sigma=diag(params[idx:(idx+n-1)])
    idx = idx + n
  
    corr = matrix(rep(0,n*n), ncol = n)
    k=1
    for(i in 1:n)
      for(j in i:n){
        if (j!=i){
          corr[i,j]=params[idx+k-1]
          corr[j,i]=params[idx+k-1]
          k=k+1
        }
        else
          corr[i,j]=1
      }
    
    idx = idx + n*(n-1)/2
    
    S = sigma %*% corr %*% sigma
    
    theta = params[idx:(idx+n-1)]
    idx = idx+n
    
    delta = params[idx:(idx+n-1)]
    idx = idx+n
    
    lambda = params[idx:(idx+n-1)]
    idx = idx+n
  }
  
  corr_matrix = matrix(ncol = n, nrow = n)
  
  corr_matrix = corr

  return(list( mu = mu, sigma =diag(sigma), corr = corr_matrix, theta = theta, delta = delta, lambda =lambda, S = S))
  
}
