# auxiliary functions

# Plots three graphs for the returns
plot_returns=function(dates, returns_daily, name){
  
  if(dates[1]>dates[2]){
    dates=rev(dates)
    returns_daily = rev(returns_daily)
  }
  
  x11()
  layout(matrix(c(1,1,2,2,1,1,2,2,3,3,2,2), 3,4, byrow = TRUE))
  plot(dates, cumsum(returns_daily), type = 'l', main=paste('log-return of',name))
  grid()
  hist(returns_daily, breaks=50, freq = FALSE, main=paste('histogram of',name))
  xx= seq(min(returns_daily), max(returns_daily), length.out = 1000)
  lines(xx, dnorm(xx, mean=mean(returns_daily), sd = sd(returns_daily)))
  legend('topright', legend="Gaussian Interp", lty = 1, col = 'black')
  plot(dates, returns_daily)
  abline(h=0, col='black')
  grid()
  
  
}


# Aux to move from array to list representation
list_from_array = function(params, model = "heston"){
  
  if(length(params)>5 & model=="heston"){
    warning("More parameters than necessary passed to the function. Only using first 5 for Heston model.")
  }
  else if(length(params)>8 & model=="bates"){
    warning("More parameters than necessary passed to the function. Only using first 8 for Bates model.")
  }
  
  r = params[1]
  k=params[2]
  eta=params[3]
  theta=params[4]
  rho = params[5]
  
  feller = k*eta*2/theta^2
  
  if (feller< 1){
    k = theta^2 / (2*eta)+ 1e-6
    warning(paste("Feller condition not verified, changing value of k to make it true. New k:",k, "New Feller value:", k*eta*2/theta^2))
  }
  
  ret_list = list(r = r, k=k, eta=eta, theta=theta, rho = rho)
  if (model=="bates"){
    ret_list["mu_j"] = params[6]
    ret_list["sigma_j"] = params[7]
    ret_list["lambda"] = params[8]
  }
  
  return(ret_list)
}


# Aux to move from array to list representation
array_from_list = function(param_list, model="heston"){
  if(!all(c('r', 'k','eta', 'theta', 'rho') %in% names(param_list))){
    stop("Missing some of the parameter in the list. Make sure it has r, k, eta, theta and rho in it.")
  }
  if (model=='bates' & !all(c('mu_j','sigma_j', 'lambda') %in% names(param_list))){
    stop("Missing some of the jump parameter in the list. Make sure it has mu_j, sigma_j and lambda in it.")
  }
  
  res = c(param_list$r, param_list$k, param_list$eta, param_list$theta, param_list$rho)
  
  
  if(model=="bates"){
    res = c(res, c(param_list$mu_j, param_list$sigma_j, param_list$lambda))
  }
  
  return(res)
}


# Function to check Feller condition and modify parameters if it is not satisfied
check_feller= function(params, eps = 1e-6){
  if (typeof(params)=="list"){
    params = array_from_list(params)
  }
  
  r = params[1]
  k=params[2]
  eta=params[3]
  theta=params[4]
  rho = params[5]
  
  feller = 2*k*eta/theta^2
  if (feller< 1){
    # rather than simply changing the value of either k or theta
    # we modify both to have a smaller difference with original values
    new_theta = sqrt(2*k*eta) 
    final_theta = (theta + new_theta)*0.5
    
    final_k = final_theta^2 /(2*eta) + eps
    
    if (2*final_k*eta <final_theta^2){
      stop("Something went terribly wrong.")
    }
    else {
      warning(paste("Feller condition not verified, changing value of k and theta to make it true. New k:",final_k, 
                    "New theta:", final_theta, "New Feller value:", final_k*eta*2/final_theta^2))
    }
    
    # saving changes to params
    if(typeof(params)=="list"){
      params$k = final_k
      params$theta = final_theta
    }
    else {
      params[2]=final_k
      params[4]=final_theta
    }
    
  }
  else{
    print(paste("Feller is verified with value", feller))
  }
  
  return(params)
}


# Functions to compute relevant moments 

mean_from_pdf = function(x, pdf, dx =x[2]-x[1]){
  return(sum(pdf*dx*x))
}


variance_from_pdf = function(x, pdf, dx =x[2]-x[1]){
  m1 = sum(pdf*dx*x)
  return(sum(pdf*dx*(x-m1)^2))
}


skewness_from_pdf = function(x, pdf, dx =x[2]-x[1]){
  
  m1 = sum(pdf*dx*x)
  m2 = sum(pdf*dx*(x-m1)^2)
  m3 = sum(pdf*dx*(x-m1)^3)
  
  skew = m3/m2**(1.5)
  return(skew)
}


kurtosis_from_pdf = function(x, pdf, dx =x[2]-x[1]){
  
  m1 = sum(pdf*dx*x)
  m2 = sum(pdf*dx*(x-m1)^2)
  m4 = sum(pdf*dx*(x-m1)^4)
  
  kurt = m4/m2**(2)
  return(kurt)
}
