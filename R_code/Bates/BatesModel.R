#####################################
##### BATES and HESTON models #######
#####################################

# based on the adaptation of some function from 
# library(NMOF)

source("ex_bates_heston.R")

library(cubature)
library(pracma)

pdfHeston = function(x,x_0, dt, sigma_0, r, k,eta, theta, rho, lower=0, upper=10, check_feller=TRUE){
  l= length(x)
  if(length(dt)!=l) stop("Time intervals and returns vectors should have same length.")
  
  y= rep(0,l)
  for(i in 1:l){
    y[i] = integrate(f = integrand_H, x = x[i], x_0 = x_0, tau = dt[i],r = r,v0 = sigma_0^2,vT = eta, rho = rho,k = k,sigma = theta,
                     lower = lower, upper=upper, abs.tol = 1e-5)$value/pi
    
    
    # y[i] = quadinf(f = integrand_H, x = x[i], x_0 = x_0, tau = dt[i],r = r,v0 = sigma_0^2,vT = eta, rho = rho,k = k,sigma = eta,
    #                  xa = lower, xb=upper, tol = 1e-5)$Q/pi
    if(y[i]<1e-8){
      y[i]=1e-8
    }
  }
  
  # # include feller condition
  # if (check_feller)
  #   if(2*k*eta<theta^2) y= rep(1e8,l)
  return(y)
}


integrand_H = function(u, x, x_0, tau ,r,v0,vT,rho,k,sigma){
  cf= my_cfHeston(om = u, x_0=x_0 ,tau = tau,r = r,v0 = sigma_0^2,vT = eta,rho = rho,k = k,sigma = theta)
  expon = exp(-1i*x*u)
  #print(paste("cf= ",cf, ", expon= ",expon))
  
  return(Re(cf*expon))
}

# Neg logLikelihood function
negloglikHeston = function(params, x, x_0, sigma_0, dt){
  pdfs= pdfHeston(x=x,dt=dt, x_0=x_0, sigma_0 = sigma_0,
                  r = params[1],k=params[2],eta=params[3],theta = params[4],rho = params[5])

  to_sum = log(pdfs)
  # print(cbind(pdfs,to_sum))
  nll = -sum(to_sum) 
    if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
      nll = 1e10
    }
  # FELLER CONDITION
  # 2*k*eta > theta^2
  # print(2*params[2]*params[3]>params[4]^2)
  return(nll)
}


#########################################################
###################### BATES ############################
#########################################################

pdfBates = function(x,x_0, dt, sigma_0, r, k,eta, theta, rho,lambda, mu_j, sigma_j, lower=0, upper=10){
  l= length(x)
  if(length(dt)!=l) stop("Time intervals and returns vectors should have same length.")
  
  y= rep(0,l)
  for(i in 1:l){
    y[i] = integrate(f = integrand_B, x = x[i], x_0 = x_0, tau = dt[i], r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
                     lambda = lambda, mu_j = mu_j, sigma_j = sigma_j,
                     lower = lower, upper=upper,abs.tol = 1e-5)$value/pi
    
    if(is.infinite(y[i])) stop(paste("Integral is infinite at ", i))
    
    if(y[i]<1e-8){
      y[i]=1e-8
    }
  }
  return(y)
}


integrand_B = function(u, x, x_0, tau ,r ,v0,vT,rho,k,sigma,lambda, mu_j, sigma_j){
  cf = my_cfBates(om = u, x_0=x_0 ,tau = tau,r = r,v0 = sigma_0^2,vT = eta, rho = rho, k = k, sigma = theta,
                  lambda = lambda, muJ = mu_j, vJ = sigma_j^2)
  expon=exp(-1i*x*u)
  
  #print(paste("cf= ",cf, ", expon= ",expon))
  return(Re(cf*expon))
}

# Neg logLikelihood function
negloglikBates = function(params, x, x_0, sigma_0, dt){
  pdfs= pdfBates(x=x,dt=dt, x_0=x_0, sigma_0 = sigma_0,
                  r = params[1],k=params[2],eta=params[3],theta = params[4],rho = params[5],
                  lambda=params[6],mu_j= params[7], sigma_j=params[8])
  
  to_sum = log(pdfs)
  # print(cbind(pdfs,to_sum))
  nll = -sum(to_sum) 
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }
  # FELLER CONDITION
  # 2*k*eta > theta^2
  # print(2*params[2]*params[3]>params[4]^2)
  return(nll)
}
