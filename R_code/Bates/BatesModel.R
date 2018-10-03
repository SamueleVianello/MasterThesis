#####################################
##### BATES and HESTON models #######
#####################################

# based on the adaptation of some function from 
# library(NMOF)


pdfHeston = function(x, dt, x_0, sigma_0, r, k,eta, theta, rho, lower=0, upper=50){
  l= length(x)
  y= rep(0,l)
  for(i in 2:l){
    y[i] = integrate(f = integrand, x = x[i], x_0 = x_0, tau = dt,r = r,v0 = sigma_0^2,vT = eta,rho = rho,k = k,sigma = eta,
                     lower = lower, upper=upper)$value/pi
    if(y[i]<1e-8){
      y[i]=1e-8
    }
  }
  return(y)
}


# Characteristic function of the Heston Model
my_cfHeston = function (om, x_0, tau, r, v0, vT, rho, k, sigma)
{
  if (sigma < 1e-08) 
    sigma <- 1e-08
  d <- sqrt((rho * sigma * (0+1i) * om - k)^2 + sigma^2 * ((0+1i) * om + om^2))
  
  g <- (k - rho * sigma * (0+1i) * om - d)/(k - rho * sigma * (0+1i) * om + d)
  
  cf1 <- (0+1i) * om * (x_0 + (r) * tau) 
  
  cf2 <- vT * k/(sigma^2) * ((k - rho * sigma * (0+1i) * om - d) * tau - 2 * log((1 - g * exp(-d * tau))/(1 - g)))
  
  cf3 <- v0/sigma^2 * (k - rho * sigma * (0+1i) * om - d) * (1 - exp(-d * tau))/(1 - g * exp(-d * tau))
  exp(cf1 + cf2 + cf3)
}

integrand = function(u, x, x_0, tau ,r , q ,v0,vT,rho,k,sigma){
  return(Re(my_cfHeston(om = u, x_0=x_0 ,tau = tau,r = r,v0 = sigma_0^2,vT = eta,rho = rho,k = k,sigma = eta)
            * exp(-1i*x*u)))
}

# Neg logLikelihood function
negloglikHeston = function(params, x, x_0, sigma_0, dt){
  pdfs= pdfHeston(x=x,dt=dt, x_0=x_0, sigma_0 = sigma_0,
                  r = params[1],k=params[2],eta=params[3],theta = params[4],rho = params[5])

  to_sum = log(pdfs)
  print(cbind(pdfs,to_sum))
  nll = -sum(to_sum) 
    if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
      nll = 1e10
    }
  return(nll)
}
