#####################################
##### BATES and HESTON models #######
#####################################

# based on the adaptation of some function from 
# library(NMOF)


library(cubature)
library(pracma)
library(gaussquad)

#########################################################
###################### HESTON ###########################
#########################################################


pdfHeston = function(x,x_0, dt, sigma_0, r, k,eta, theta, rho, lower=0, upper=50, check_feller=TRUE){
  l= length(x)
  if(length(dt)!=l) stop("Time intervals and returns vectors should have same length.")
  
  y= rep(0,l)
  for(i in 1:l){
    # y[i] = integrate(f = integrand_H, x = x[i], x_0 = x_0, tau = dt[i],r = r,v0 = sigma_0^2,vT = eta, rho = rho,k = k,sigma = theta,
    #                  lower = lower, upper=upper, abs.tol = 1e-5)$value/pi
    
    # y[i] = gauss_kronrod(f = integrand_H, x = x[i], x_0 = x_0, tau = dt[i],r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
    #                  a =lower, b=upper)$value/pi ####### WRONG RESULTS!
    
    # y[i] = quadinf(f = integrand_H, x = x[i], x_0 = x_0, tau = dt[i],r = r,v0 = sigma_0^2,vT = eta, rho = rho,k = k,sigma = eta,
    #                  xa = lower, xb=upper, tol = 1e-5)$Q/pi ###### TOO LONG
    
    y[i] = clenshaw_curtis(f = integrand_H, x = x[i], x_0 = x_0, tau = dt[i],r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
                         a =lower, b=upper, n=2**8)/pi
    
    
    
    # print(y[i])
    if(is.nan(y[i])|| y[i]<1e-8 ){
      y[i]=1e-8
    }
  }
  
  # # include feller condition
  # if (check_feller)
  #   if(2*k*eta<theta^2) y= rep(1e8,l)
  return(y)
}


my_cfHeston= function (om, x_0, tau, r, v0, vT, rho, k, sigma) 
{
  if (sigma < 1e-08) 
    sigma <- 1e-08
  
  om[which(om==0)]= 1e-8
  d <- sqrt((rho * sigma * (0+1i) * om - k)^2 + sigma^2 * ((0+1i) *om + om^2))
  g <- (k - rho * sigma * (0+1i) * om - d)/(k - rho * sigma * (0+1i) * om + d)
  
  # to avoid nan due to 0/0
  g[which((k - rho * sigma * (0+1i) * om - d)==0)]=0
  
  cf1 <- (0+1i) * om * (x_0 + r * tau)
  
  cf2 <- vT * k/(sigma^2) * ((k - rho * sigma * (0+1i) * om - d) * tau - 2 * log((1 - g * exp(-d * tau))/(1 - g)))
  
  cf3 <- v0/sigma^2 * (k - rho * sigma * (0+1i) * om - d) * (1 - exp(-d * tau))/(1 - g * exp(-d * tau))
  
  # print(paste("exp1=", cf1, ", exp2=",cf2, ", exp3=",cf3))
  res = exp(cf1 + cf2 + cf3)
  if (is.na(res) || is.nan(res)){
    print(paste(om,x_0,tau,r,v0,vT,rho,k,sigma))
    print(paste("exp1=", cf1, ", exp2=",cf2, ", exp3=",cf3, "d=",d, "g=",g))
  }
  return(exp(cf1 + cf2 + cf3))
}


integrand_H = function(u, x, x_0, tau ,r,v0,vT,rho,k,sigma){
  cf= my_cfHeston(om = u, x_0=x_0 ,tau = tau,r = r,v0 = v0,vT = vT,rho = rho,k = k,sigma = sigma)
  expon = exp(-1i*x*u)
  # print(paste("cf= ",cf, ", expon= ",expon))
  # print(paste("HI!",Re(cf*expon)))
  return(Re(cf*expon))
}

# Neg logLikelihood function
negloglikHeston = function(params, x, x_0, sigma_0, dt, model="heston_ab"){
  if (model=='heston'){
    # using k and eta
    pdfs= pdfHeston(x=x, dt=dt, x_0=x_0, sigma_0 = sigma_0,
                    r = params[1],k=params[2], eta=params[3], theta = params[4], rho = params[5])
  }
  else if(model=="heston_ab"){
    # using alpha and beta
    pdfs= pdfHeston(x=x,dt=dt, x_0=x_0, sigma_0 = sigma_0,
                    r = params[1],k=params[3], eta=params[2]/params[3], theta = params[4], rho = params[5])
  }
  

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

loglikHeston = function(params, x, x_0, sigma_0, dt){
  pdfs= pdfHeston(x=x,dt=dt, x_0=x_0, sigma_0 = sigma_0,
                  r = params[1],k=params[2], eta=params[3], theta = params[4], rho = params[5])
  
  to_sum = log(pdfs)
  # print(cbind(pdfs,to_sum))
  nll = sum(to_sum) 
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

pdfBates = function(x, x_0, dt, sigma_0, r, k,eta, theta, rho,lambda, mu_j, sigma_j, lower=0, upper=50){
  l= length(x)
  if(length(dt)!=l) stop("Time intervals and returns vectors should have same length.")
  
  y= rep(0,l)
  for(i in 1:l){
    # y[i] = integrate(f = integrand_B, x = x[i], x_0 = x_0, tau = dt[i], r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
    #                  lambda = lambda, mu_j = mu_j, sigma_j = sigma_j,
    #                  lower = lower, upper=upper,abs.tol = 1e-6)$value/pi
    
    # print("Using gauss kronrod.")
    # y[i] = gauss_kronrod(f = integrand_B, x = x[i], x_0 = x_0, tau = dt[i], r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
    #                  lambda = lambda, mu_j = mu_j, sigma_j = sigma_j,
    #                  a = lower, b=upper)$value/pi

    y[i] = clenshaw_curtis(f = integrand_B, x = x[i], x_0 = x_0, tau = dt[i],r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
                           lambda = lambda, mu_j = mu_j, sigma_j = sigma_j,
                           a =lower, b=upper, n=2**8)/pi
    
    # rulez =glaguerre.quadrature.rules(30, alpha =0.5)[[30]]
    # y[i]= glaguerre.quadrature(functn = integrand_B,alpha =0.5, rule = rulez,x = x[i], x_0 = x_0, tau = dt[i],r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
    #                              lambda = lambda, mu_j = mu_j, sigma_j = sigma_j,
    #                              lower = lower,upper = Inf)/pi

    if(is.infinite(y[i])) stop(paste("Integral is infinite at ", i))
    
    if(y[i]<1e-8){
      y[i]=1e-8
    }
  }
  return(y)
}


integrand_B = function(u, x, x_0, tau ,r ,v0,vT,rho,k,sigma,lambda, mu_j, sigma_j){
  cf = my_cfBates(om = u, x_0=x_0 ,tau = tau,r = r,v0 = v0,vT = vT, rho = rho, k = k, sigma = sigma,
                  lambda = lambda, muJ = mu_j, vJ = sigma_j^2)
  expon=exp(-1i*x*u)
  
  #print(paste("cf= ",cf, ", expon= ",expon))
  return(Re(cf*expon))
}

# Neg logLikelihood function
negloglikBates = function(params, x, x_0, sigma_0, dt, model="bates_ab"){
  if (model=="bates"){
    # using k and eta
    pdfs= pdfBates(x=x,dt=dt, x_0=x_0, sigma_0 = sigma_0,
                   r = params[1],k=params[2],eta=params[3],theta = params[4],rho = params[5],
                   lambda=params[6],mu_j= params[7], sigma_j=params[8])
  }
  else if(model=="bates_ab"){
    # using alpha and beta
    pdfs= pdfBates(x=x,dt=dt, x_0=x_0, sigma_0 = sigma_0,
                   r = params[1],k=params[3],eta=params[2]/params[3],theta = params[4],rho = params[5],
                   lambda=params[6],mu_j= params[7], sigma_j=params[8])
  }

  
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




my_cfBates= function (om, x_0, tau, r, v0, vT, rho, k, sigma, lambda, muJ, vJ) 
{
  if (sigma < 1e-08) 
    sigma <- 1e-08
  #sigma <- max(sigma, 1e-04)
  
  om[which(om==0)]= 1e-9
  om1i <- om * (0+1i)
  
  d <- sqrt((rho * sigma * om1i - k)^2 + sigma^2 * (om1i + om^2))
  
  g <- (k - rho * sigma * om1i - d)/(k - rho * sigma * om1i + d)
  
  # to avoid nan due to 0/0
  g[which((k - rho * sigma * (0+1i) * om - d)==0)]=0
  
  cf1 <- om1i * (x_0 + r * tau) #ok
  
  cf2 <- vT * k/(sigma^2) * ((k - rho * sigma * om1i - d) * tau - 2 * log((1 - g * exp(-d * tau))/(1 - g)))
  
  cf3 <- v0/sigma^2 * (k - rho * sigma * om1i - d) * (1 - exp(-d * tau))/(1 - g * exp(-d * tau))
  
  cf4 <- -lambda * muJ * om1i * tau + lambda * tau * ((1 + muJ)^(om1i) * exp(vJ * (om1i/2) * (om1i - 1)) - 1)
  
  
  res= exp(cf1 + cf2 + cf3 + cf4)
  
  if (is.na(res) || is.nan(res)){
    print(paste(om,x_0,tau,r,v0,vT,rho,k,sigma, lambda, muJ,vJ))
    print(paste("exp1=", cf1, ", exp2=",cf2, ", exp3=",cf3,"exp4=",cf4, "d=",d, "g=",g))
  }
  return(res)
}
