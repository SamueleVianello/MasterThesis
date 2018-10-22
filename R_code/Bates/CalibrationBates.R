###########################################
##### Calibration of Bates and Heston #####
###########################################

# need cumulative returns X_t =  log (S_t/S_0)

CalibrateModel=function(x, x_0,sigma_0, dt, trace = 10, initial, deoptim=FALSE, model="heston_ab", 
                        feller = TRUE, sigma_is_param = TRUE){
  
  library(DEoptim)
  library(GenSA)
  source("BatesModel.R")
  
  bounds = BoundsCreator(model = model, sigma_param = sigma_is_param)
  if(model =='heston_ab'||model=="heston") {obj=negloglikHeston}
  else if (model =='bates') {obj=negloglikBates}
  else{ stop("Choose model between heston or bates")}
  
  if (sigma_is_param){
    # First optimization using deoptim
    if(deoptim){
      start_time_deoptim <- Sys.time()
      # print("Starting calibration using DEoptim...")
      # control_list_deoptim = list(itermax = 10, NP = 200, strategy = 6,trace=trace)
      # 
      # 
      # outDE <- DEoptim(obj,
      #                  lower = bounds$lower, upper = bounds$upper, control = control_list_deoptim,
      #                  x=x, x_0 = x_0, dt = dt, model=model, check_feller = feller)
      # 
      # 
      # initial=outDE$optim$bestmem
      
      
      
      outSA= GenSA(fn = negloglikHeston, lower = bounds$lower, upper = bounds$upper, par = initial,
                             control = c(max.time=300,verbose=TRUE, simple.function=FALSE),
                             x = x, x_0 = x_0, dt = dt, model= model)
      
      initial=outSA$par
      
      end_time_deoptim <- Sys.time()
      
    }
    # Second and final optimization using nlminb
    print("Starting calibration using nlminb...")
    
    
    start_time_nlminb <- Sys.time()
    out_nlminb = nlminb(initial,objective = obj,lower = bounds$lower, upper = bounds$upper,
                        x=x, x_0 = x_0,dt = dt, model=model,  check_feller = feller,
                        control=list(eval.max = 1000,iter.max = 100, trace = trace))
    print(out_nlminb)
    end_time_nlminb <- Sys.time()
  }
  
  
  else{ # sigma is not a parameter
    # First optimization using deoptim
    if(deoptim){
      start_time_deoptim <- Sys.time()
      # print("Starting calibration using DEoptim...")
      # control_list_deoptim = list(itermax = 10, NP = 200, strategy = 6,trace=trace)
      # 
      # 
      # outDE <- DEoptim(obj,
      #                  lower = bounds$lower, upper = bounds$upper, control = control_list_deoptim,
      #                  x=x, x_0 = x_0, sigma_0=sigma_0, dt = dt, model=model, check_feller = feller)
      # 
      # 
      # initial=outDE$optim$bestmem
      
      outSA= GenSA(fn = negloglikHeston, lower = bounds$lower, upper = bounds$upper, par = initial,
                   control = c(max.time=300,verbose=TRUE,simple.function=TRUE),
                   x = x, x_0 = x_0, dt = dt,sigma_0=sigma_0, model= model)
      
      initial=outSA$par
      
      end_time_deoptim <- Sys.time()
    }
    # Second and final optimization using nlminb
    print("Starting calibration using nlminb...")
    
    
    start_time_nlminb <- Sys.time()
    out_nlminb = nlminb(initial,objective = obj,lower = bounds$lower, upper = bounds$upper,
                        x=x, x_0 = x_0, sigma_0=sigma_0,dt = dt, model=model,  check_feller = feller,
                        control=list(eval.max = 1000,iter.max = 100, trace = trace))
    print(out_nlminb)
    end_time_nlminb <- Sys.time()
  }
  
  if(deoptim){
    print(paste("DEoptim:",end_time_deoptim - start_time_deoptim))
  }
  print(paste("nlminb:", end_time_nlminb - start_time_nlminb))
  
  if(deoptim){
    print(paste("TOTAL:", end_time_nlminb - start_time_deoptim))
  }
  
  res = ParametersReconstruction(out_nlminb$par, model = model, sigma_param = sigma_is_param)
  res[["message"]] = out_nlminb$message
  res[["objective_function"]] = out_nlminb$objective
  res[["total_time"]] = end_time_nlminb - start_time_nlminb
  return(res)
}






BoundsCreator= function(n=1, model = "heston_ab", sigma_param){
  # Creates lower and upper boundaries for the DEoptim optimization on the likelihood
  # for a n-multivariate merton process and n_common common jumps
  min_mu = -3
  max_mu = 3
  min_lambda = 0.001
  max_lambda = 100
  min_sigma = 1e-5
  max_sigma = 2
  min_corr = -1
  max_corr = 1

  if (model=="heston"|| model=="heston_ab"){
    # r, a, b, theta, rho
    low = c(min_mu,0.00001,min_sigma,min_sigma,min_corr)
    up = c(max_mu,max_mu,max_sigma,1,max_corr)
    
    if (sigma_param){
      low= c(low,min_sigma)
      up = c(up, max_sigma)
    }
  }
  else if(model=="bates"||model=="bates_ab") {
    # r, k, eta, theta, rho, lambda, mu_j, sigma_j 
    low = c(min_mu,0.00001,min_sigma,min_sigma,min_corr,min_lambda, min_mu, min_sigma)
    up = c(max_mu,max_mu,max_sigma,1,max_corr, max_lambda, max_mu, max_sigma)
  }
  

  return(list(lower = low, upper = up))
}





ParametersReconstruction = function(params, model="heston", sigma_param){
  if(model =="heston"){
    print("Model is heston")
    res=list(r=params[1],k=params[2],eta=params[3],theta=params[4],rho=params[5])

  }
  else if (model =="heston_ab"){
    print("Model is heston_ab")
    res=list(r=params[1],k=params[3],eta=params[2]/params[3],theta=params[4],rho=params[5])
  }
  else if (model == "bates"){
    res=list(r=params[1],k=params[2],eta=params[3],theta=params[4],rho=params[5], lambda = params[6], mu_j=params[7],sigma_j=params[8])
  }
  else if (model == "bates_ab"){
    res=list(r=params[1],k=params[3],eta=params[2]/params[3],theta=params[4],rho=params[5], lambda = params[6], mu_j=params[7],sigma_j=params[8])
  }
  
  if (sigma_param){
    res["sigma_0"]=params[6]
  }
  return(res)
}






