###########################################
##### Calibration of Bates and Heston #####
#####   based on maximum likelihood   #####
###########################################

# library(adagio)
# library(GenSA)
library(DEoptimR)

source("BatesModel.R")
source("AuxiliaryFunctions.R")


CalibrateModel=function(x, dt,   # log-returns from data and the time interval length in years
                        initial, # initial guess for parameters
                        model,   # 'heston' or 'bates'
                        global=TRUE,local=TRUE, # whether to perform global and local optimization
                        feller = TRUE,          # if we want to make sure the feller condition is verified
                        trace = 10,             # controls the print output of optimizers
                        fixing=NA,                         # (used for drift here) which parameters to keep fixed
                        fixing_values = initial,           # (used for drift here) to which value we should fix the fixed parameters
                        set_upper = NA, upper_bnd = NA,    # see BoundsCreator function
                        set_lower = NA, lower_bnd=NA){     # see BoundsCreator function 
  
  
  
  print(paste("Initial: ", paste(initial, collapse = "    ")))
  #print(fixing)
  bounds = BoundsCreator(model=model, fixing = fixing, fixing_values = fixing_values,
                         set_upper =set_upper, upper_bnd=upper_bnd,
                         set_lower = set_lower, lower_bnd = lower_bnd)
  print(bounds)
  
  control_GenSA = c(maxit=1000,max.time=60, verbose=TRUE, simple.function=FALSE)
  #extra param: maxit=1000, max.time=60, 
  control_nlminb = list(eval.max = 10000,iter.max = 500, trace = trace, rel.tol = 1e-10, abs.tol=0, x.tol=1e-10)
  
  if(model=="heston"){
    obj = negloglikHeston
  }
  else if (model=="bates"){
    obj = negloglikBates
  }
  else {
    stop("Please choose model between heston or bates.")
  }
  
  
  if(global){
    # Optimization using global optimizer
    
    # NOTE: Using the negloglik as objective function, the best optimizer was
    # JDEoptim since it allows for non linear constraints (feller condition)
    
    # print("Starting calibration using GenSA...") # make sure to add the GenSA library
    # outSA= GenSA(fn = negloglikBates, lower = bounds$lower, upper = bounds$upper, par = initial,
    #              control = control_GenSA,
    #              x = x, dt = dt, check_feller = feller)
    # 
    # initial=outSA$par
    # res_param = outSA$xmin
    # res_obj = outSA$value
    
    print("Starting calibration using JDEoptim...")
    
    feller_constr = function(params,... ){
        k=params[2]
        eta=params[3]
        theta=params[4]
        
        feller = 2*k*eta/theta^2 -1
        
        return(-feller) # because inequality has to be of the form f(x)<=0 
      }
    start_time_global <- Sys.time()
    outJDE = JDEoptim(fn = obj, lower = bounds$lower, upper = bounds$upper,
                      constr=feller_constr, meq=0,
                      trace = trace, triter=trace, tol = 1, maxiter = 1000, NP = 20*length(bounds$lower),
                      x=x, dt = dt, check_feller = F)
                      # try with tol = 1
    end_time_global <- Sys.time()
    print(outJDE)
    initial = outJDE$par
    res_param= outJDE$par
    res_obj=outJDE$value
    
    
    # print("Starting calibration using Nelder-Mead...") # make sure to add the adagio library
    # start_time_global <- Sys.time()
    # outNM = nelminb(fn = obj, x0 = initial, lower = bounds$lower, upper = bounds$upper,
    #                 tol = 1e-6, step = rep(0.0001, length(bounds$lower)),maxfeval = 1e5,
    #                 x = x, dt = dt, check_feller = feller)
    # end_time_global <- Sys.time()
    # print(outNM)
    # initial = outNM$xmin
    # res_param = outNM$xmin
    # res_obj = outNM$fmin
    
  }
  
  if(local){
    # Local optimization using nlminb
    print("Starting calibration using nlminb...")
    
    
    start_time_local <- Sys.time()
    out_nlminb = nlminb(initial,objective = obj, lower = bounds$lower, upper = bounds$upper,
                        x=x, dt = dt,  check_feller = feller,
                        control=control_nlminb)
    print(out_nlminb)
    
    res_param = out_nlminb$par
    res_obj = out_nlminb$objective
    end_time_local <- Sys.time()
    
  }
  
  if(global){
    print(paste("Global optimization time:",end_time_global - start_time_global))
  }
  
  if(local){
    print(paste("Local optimization time:", end_time_local - start_time_local))
  }
  
  if(global & local){
    print(paste("TOTAL optimization time:", end_time_local - start_time_global))
  }
  
  res = ParametersReconstruction(res_param, model=model)
  res[["objective_function"]] = res_obj
  # res[["total_time"]] = end_time_local - ifelse(global, start_time_global,start_time_local)
  # res[["message"]] = out_nlminb$message
  return(res)
}





BoundsCreator= function(model = "heston",   # 'heston' or 'bates'
                        fixing=NA,          # which parameters to fix
                        fixing_values=NA,   # value to which fix range 
                        eps=0.05,           # free range for fixed params
                        set_upper = NA,     # which upper bounds to custom set
                        upper_bnd = NA,     # values of custom upper bounds
                        set_lower = NA,     # which lower bounds to custom set
                        lower_bnd=NA){      # values of custom lower bounds
  
  # Creates lower and upper boundaries for the calibration paramenters
  
  # checking inputs
  if( (is.na(set_upper[1]) & !is.na(upper_bnd[1])) | (!is.na(set_upper[1]) & is.na(upper_bnd[1])) | length(set_upper)!=length(upper_bnd) ){
    stop("Input values for upper boundary are not consistent.")
  }
  if( (is.na(set_lower[1]) & !is.na(lower_bnd[1])) | (!is.na(set_lower[1]) & is.na(lower_bnd[1])) | length(set_lower)!=length(lower_bnd) ){
    stop("Input values for lower boundary are not consistent.")
  }
  if((!is.na(fixing[1]) & is.na(fixing_values[1])) ){
    stop("Input for values to be fixed are not consistent.")
  }
  
  min_mu = -3
  max_mu = 3
  
  min_k = 1e-5
  max_k = 2
  min_eta = 1e-5
  max_eta = 3
  min_theta = 1e-5
  max_theta = 2
  
  min_corr = -1
  max_corr = -0.0001

  # min_muj = -0.5
  max_muj = -0.2541196
  min_muj =-0.5
  #max_muj = -0.028
  min_sigmaj = 1e-5
  max_sigmaj = 0.1
  min_lambda = 1e-5
  max_lambda = 10

  # beta = k 
  print(paste("Model is", model))
  
  if (model=="heston"){
    # r, a, b, theta, rho
    low = c(min_mu,min_k,min_eta,min_theta,min_corr)
    up = c(max_mu,max_k,max_eta,max_theta,max_corr)
  }
  else if(model=="bates") {
    # r, k, eta, theta, rho,  mu_j, sigma_j, lambda
    low = c(min_mu,min_k,min_eta,min_theta,min_corr, min_muj, min_sigmaj,min_lambda)
    up = c(max_mu,max_k,max_eta,max_theta,max_corr, max_muj, max_sigmaj, max_lambda)
  }
  
  if ( !is.na(fixing[1]) & !is.na(fixing_values[1])){
    print(paste("We are fixing parameters: ", paste(fixing,collapse=","), "with free range of +/-", eps))
    for (i in 1:length(fixing)){
      low[fixing[i]]=fixing_values[fixing[i]] - eps
      up[fixing[i]]=fixing_values[fixing[i]]+ eps
    }
  }
  else{
    print("No parameters are being fixed")
  }
  
  
  if( !is.na(set_upper[1]) & !is.na(upper_bnd[1])){
    for(i in 1:length(set_upper)){
      print(set_upper[i])
      up[set_upper[i]]=upper_bnd[i]
    }
  }
  
  if( !is.na(set_lower[1]) & !is.na(lower_bnd[1])){
    for(i in 1:length(set_lower)){
      low[set_lower[i]]=lower_bnd[i]
    }
  }
  
  
  return(list(lower = low, upper = up))
}





ParametersReconstruction = function(params, model="heston"){
  if(model =="heston"){
    print("Model is heston")
    res=list(r=params[1],k=params[2],eta=params[3],theta=params[4],rho=params[5])

  }
  else if (model == "bates"){
    res=list(r=params[1],k=params[2],eta=params[3],theta=params[4],rho=params[5], mu_j=params[6],sigma_j=params[7], lambda = params[8])
  }
  return(res)
}
