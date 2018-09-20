###############################################################
############ Calibration of a Multivariate Merton #############
###############################################################


CalibrateMVMerton=function(x, n, dt, trace = 10){
  
  library(DEoptim)
  source("MultivariateMertonModel.R")
  
  # Choose objective function given the number of assets
  if(n ==1) obj = negloglik_1asset_nocommon
  else if(n ==2) obj = negloglik_2assets_nocommon
  else if(n ==3) obj = negloglik_3assets_nocommon
  else if(n ==4) obj = negloglik_4assets_nocommon
  else obj = negloglik
  
  bounds_nocommon = BoundsCreator(n, n_common=0)
  
  # First optimization using deoptim
  print("Starting calibration using DEoptim...")
  control_list_deoptim = list(itermax = 50, NP = 200, strategy = 6,trace=trace)
  
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
  
  
  return(ParametersReconstruction(out_nlminb$par,n=n,common = FALSE))
}






BoundsCreator= function(n, n_common=1 ){
  # Creates lower and upper boundaries for the DEoptim optimization on the likelihood
  # for a n-multivariate merton process and n_common common jumps
  min_mu = -10
  max_mu = 10
  min_lambda = 0.1
  max_lambda = 100
  min_theta = -1
  max_theta = -0.1 
  min_var = 1e-4
  max_var = 10
  min_alpha = -1
  max_alpha = 1
  
  # initialising resulting vector low and up
  leng = 4*n + (n+1)*n*0.5 + n_common*n + 3*n_common
  low = rep(0,leng)
  up = rep(0,leng)
  
  
  idx =1
  # mean of continuos part
  low[idx:(idx+n-1)] = rep(min_mu,n)
  up[idx:(idx+n-1)] = rep(max_mu,n)
  idx = idx+n 
  
  # covariance matrix of continuous part
  N_var = n*(n+1)/2
  low[idx:(idx + N_var -1)] = rep(min_var,N_var)
  up[idx:(idx + N_var -1)] = rep(max_var,N_var)
  idx = idx + N_var
  
  # means of idiosyncratic term
  low[idx:(idx+n-1)] = rep(min_theta,n)
  up[idx:(idx+n-1)] = rep(max_theta,n)
  idx = idx+n
  
  # variance of idyosincratic term
  low[idx:(idx+n-1)] = rep(min_var,n)
  up[idx:(idx+n-1)] = rep(max_var,n)
  idx = idx+n
  
  # lambda of idyosincratic poissons
  low[idx:(idx+n-1)] = rep(min_lambda,n)
  up[idx:(idx+n-1)] = rep(max_lambda,n)
  idx = idx+n
  
  if (n_common>0)
  {  
    # boundaries on the parameters of common jumps
    low[idx] =min_mu
    up[idx] = max_mu
    idx = idx+1
    
    low[idx] =min_var
    up[idx] = max_var
    idx = idx+1
    
    low[idx] =min_lambda
    up[idx] = max_lambda
    idx = idx+1
    
    # boundaries on alpha
    low[idx:(idx+n-1)] = rep(min_alpha,n)
    up[idx:(idx+n-1)] = rep(max_alpha,n)
    idx = idx+n
  }
  
  if(  ((idx-1)!=length(low))  || (length(low)!=length(up)) )
    stop("Error in parameter reconstruction: number of parameters is wrong.")
  
  return(list(lower = low, upper = up))
}


ParametersReconstruction = function(params, n, common = TRUE){
  
  # reconstruction of parameters:
  idx =1
  mu = params[idx:(idx+n-1)]
  idx = idx+n 
  
  
  S = matrix(rep(0,n*n), ncol = n)
  i=1
  j=1
  for(k in 1:(n*(n+1)/2)){
    S[i,j] = params[idx+k-1]
    S[j,i] =  S[i,j]
    j=j+1
    if(j == n+1){
      i=i+1
      j=i
    }
  }
  idx = idx + n*(n+1)/2
  
  theta = params[idx:(idx+n-1)]
  idx = idx+n
  
  delta = params[idx:(idx+n-1)]
  idx = idx+n
  
  lambda = params[idx:(idx+n-1)]
  idx = idx+n
  
  if (common){  
    theta_z = params[idx]
    idx = idx+1
    
    delta_z = params[idx]
    idx = idx+1
    
    lambda_z = params[idx]
    idx = idx+1
    
    alpha = params[idx:(idx+n-1)]
    idx = idx+n
    
    return(list( mu = mu, S = S, theta = theta, delta = delta, lambda =lambda,
                 theta_z = theta_z, delta_z = delta_z, lambda_z = lambda_z, alpha = alpha))
  }
  
  else{
    return(list( mu = mu, S = S, theta = theta, delta = delta, lambda =lambda))
  }
}
