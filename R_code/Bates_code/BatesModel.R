#####################################
##### BATES and HESTON models #######
#####################################

# based on the adaptation of some function from 
# library(NMOF)



library(pracma)


#########################################################
###################### HESTON ###########################
#########################################################

pdfHeston_fft = function(x, dt, r, k, eta, theta, rho, N = 2^10){
  
  aux_chf = function(u,tau = dt, r_cf = r, k_cf = k,vT = eta, sigma = theta, rho_cf = rho ){
    uncond_cfHeston(u,tau ,r_cf ,k_cf,vT,sigma,rho_cf)
  }
  
  min_x = min(x, -1)
  max_x = max(x, 1)
  eps = 0.01
  pdf = characteristic_function_to_density(aux_chf,N, min_x-eps, max_x+eps)
  
  pdf_xx = interp1(pdf$x, pdf$density, x,method = 'spline')
  return(pdf_xx)
}


pdfHeston = function(x, dt, r, k, eta, theta, rho, lower=0, upper=500, N=2**12){
  l= length(x)
  if(length(dt)==1){
    y = pdfHeston_fft(x,dt,r,k,eta,theta,rho,N)
  }
  else if(length(dt)==l){
    y= rep(0,l)
    for(i in 1:l){
      y[i] = clenshaw_curtis(f = integrand_uncond_H, x = x[i], tau = dt[i],r = r, vT = eta, rho = rho, k = k, sigma = theta,
                             a =lower, b=upper*5, n=N)/pi
    }
    if(is.infinite(y[i])) stop(paste("Integral is infinite at x[", i, "] = ", x[i]))
    
  }
  else {
    stop("Time intervals and returns vectors should have same length.")
  }
  
  y[which(y<1e-10)] = 1e-10
  
  return(y)
}



# unconditional characteristic function
uncond_cfHeston = function(om, tau, r, k, vT, sigma, rho){
  
  if (sigma < 1e-08) 
    sigma <- 1e-08
  
  om[which(om==0)]= 1e-8
  d <- sqrt((rho * sigma * (0+1i) * om - k)^2 + sigma^2 * ((0+1i) *om + om^2))
  g <- (k - rho * sigma * (0+1i) * om - d)/(k - rho * sigma * (0+1i) * om + d)
  
  # to avoid nan due to 0/0
  g[which((k - rho * sigma * (0+1i) * om - d)==0)]=0
  
  A <- 1i*om *r * tau + vT * k/(sigma^2) * ((k - rho * sigma * (0+1i) * om - d) * tau - 2 * log((1 - g * exp(-d * tau))/(1 - g)))
  

  
  B <- (k - rho * sigma * (0+1i) * om - d)/sigma^2 * (1 - exp(-d * tau))/(1 - g * exp(-d * tau))
  
  w = 2*k/sigma^2
  nu = 2*k*vT/sigma^2
  
  res = exp(A)*(w/(w-B))^nu
  return(res)
}

integrand_uncond_H = function(u, x, tau ,r,vT,rho,k,sigma){
  cf= uncond_cfHeston(om = u,tau = tau,r = r,  vT = vT, rho = rho, k = k, sigma = sigma)
  expon = exp(-1i*x*u)
  # print(paste("cf= ",cf, ", expon= ",expon))
  # print(paste("HI!",Re(cf*expon)))
  return(Re(cf*expon))
}


  
  
# Neg logLikelihood function
negloglikHeston = function(params, x, dt, check_feller=TRUE){
  
  r = params[1]
  k =params[2]
  eta = params[3]
  theta = params[4]
  rho = params[5]
  
  # just an extra check because NA would randomly appear from nowhere 
  if (is.na(2*k*eta*theta)){
    print(paste('check_feller=',check_feller, "2*k*eta<theta^2 ", 2*k*eta<theta^2))
    
    print(paste('k=',k,'eta=',eta, 'theta=',theta))
  }
  
  # Feller condition check
  if (check_feller && 2*k*eta<theta^2){
      nll= 1e8      # if Feller condition not verified, assign to the nll a high number
  }
  # if dont check feller or feller satisfied
  else{
    pdfs= pdfHeston(x=x, dt=dt,
                    r = r, k=k, eta=eta, theta = theta, rho = rho)
    
    to_sum = log(pdfs)
    # print(cbind(pdfs,to_sum))
    nll = -sum(to_sum) 
  }
  
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }

  return(nll)
}


# #########################################################
# ###################### BATES ############################
# #########################################################

pdfBates_fft = function(x,  dt,  r, k,eta, theta, rho, lambda, mu_j, sigma_j,N=2^12){
  
  aux_bates = function( u, tau=dt, r_cf=r, vT=eta, rho_cf=rho, k_cf=k, sigma=theta, lambda_cf=lambda, muJ=mu_j, vJ=sigma_j){
    uncond_cfBates(u,tau, r_cf, vT, rho_cf, k_cf, sigma, lambda_cf, muJ, vJ)
  }
  
  min_x = min(x, -1)
  max_x = max(x, 1)
  eps = 0.1
  pdf = characteristic_function_to_density(aux_bates,N, min_x-eps, max_x+eps)
  
  pdf_xx = interp1( x = pdf$x, y = pdf$density, xi = x, method = 'spline')
  return(pdf_xx)
}


pdfBates = function(x,  dt,  r, k,eta, theta, rho, lambda, mu_j, sigma_j, lower=0, upper=5000, N=2^12){
  l= length(x)
  if(length(dt)==1) {
    # only if all x are sample from the same dt, for instance with x = daily returns, dt = 1/255
    y = pdfBates_fft(x,  dt,  r, k, eta, theta, rho, lambda, mu_j, sigma_j,N)
  }
  else if (length(dt)==l){
    y= rep(0,l)
    for(i in 1:l){
      # y[i] = integrate(f = integrand_B, x = x[i], x_0 = x_0, tau = dt[i], r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
      #                  lambda = lambda, mu_j = mu_j, sigma_j = sigma_j,
      #                  lower = lower, upper=upper,abs.tol = 1e-6)$value/pi
      
      # print("Using gauss kronrod.")
      # y[i] = gauss_kronrod(f = integrand_B, x = x[i], x_0 = x_0, tau = dt[i], r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
      #                  lambda = lambda, mu_j = mu_j, sigma_j = sigma_j,
      #                  a = lower, b=upper)$value/pi
      
      y[i] = clenshaw_curtis(f = integrand_uncond_B, x = x[i], tau = dt[i],r = r, vT = eta, rho = rho, k = k, sigma = theta,
                             lambda = lambda, mu_j = mu_j, sigma_j = sigma_j,
                             a =lower, b=upper, n=N)/pi
      
      # rulez =glaguerre.quadrature.rules(30, alpha =0.5)[[30]]
      # y[i]= glaguerre.quadrature(functn = integrand_B,alpha =0.5, rule = rulez,x = x[i], x_0 = x_0, tau = dt[i],r = r, v0 = sigma_0^2, vT = eta, rho = rho, k = k, sigma = theta,
      #                              lambda = lambda, mu_j = mu_j, sigma_j = sigma_j,
      #                              lower = lower,upper = Inf)/pi
      
      if(is.infinite(y[i])) stop(paste("Integral is infinite at ", i))
      
    }
  }
  else {
    stop("Time intervals length are not consistent. It should be of length 1 or same length as returns.")
  }
  y[which(y<1e-10)] = 1e-10
  return(y)
}



integrand_B = function(u, x, x_0, tau ,r ,v0,vT,rho,k,sigma,lambda, mu_j, sigma_j){
  cf = my_cfBates(om = u, x_0=x_0 ,tau = tau,r = r,v0 = v0,vT = vT, rho = rho, k = k, sigma = sigma,
                  lambda = lambda, muJ = mu_j, vJ = sigma_j^2)
  expon=exp(-1i*x*u)

  #print(paste("cf= ",cf, ", expon= ",expon))
  return(Re(cf*expon))
}


integrand_uncond_B = function(u, x, tau ,r, vT, rho, k, sigma, lambda, mu_j, sigma_j){
  cf = uncond_cfBates(om = u, tau = tau,r = r, vT = vT, rho = rho, k = k, sigma = sigma,
                  lambda = lambda, muJ = mu_j, vJ = sigma_j^2)
  expon=exp(-1i*x*u)
  #print(paste("cf= ",cf, ", expon= ",expon))
  return(Re(cf*expon))
}


# Neg logLikelihood function 
negloglikBates = function(params, x, dt, model="bates", check_feller){
  r = params[1]
  k=params[2]
  eta=params[3]
  theta = params[4]
  rho = params[5]
  mu_j= params[6]
  sigma_j=params[7]
  lambda=params[8]

  if (check_feller & 2*k*eta<theta^2){
    # Check feller
    nll= 1e8
  }
  else{
    pdfs= pdfBates(x=x,dt=dt, r=r,k=k, eta=eta,theta=theta,rho =rho,
                          mu_j = mu_j, sigma_j = sigma_j, lambda=lambda)

    to_sum = log(pdfs)
    nll = -sum(to_sum)
  }
  # print(cbind(pdfs,to_sum))
  
  
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }

  return(nll)
}



# conditional on initial values of logreturns x_0 and variance V0 
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
  
  B = 1/sigma^2 * (k - rho * sigma * om1i - d) * (1 - exp(-d * tau))/(1 - g * exp(-d * tau))

  cf3 <- v0* B

  cf4 <- -lambda * muJ * om1i * tau + lambda * tau * ((1 + muJ)^(om1i) * exp(vJ * (om1i/2) * (om1i - 1)) - 1)


  res= exp(cf1 + cf2 + cf3 + cf4)

  if (is.na(res) || is.nan(res)){
    print(paste(om,x_0,tau,r,v0,vT,rho,k,sigma, lambda, muJ,vJ))
    print(paste("exp1=", cf1, ", exp2=",cf2, ", exp3=",cf3,"exp4=",cf4, "d=",d, "g=",g))
  }
  return(res)
}


uncond_cfBates= function (om,tau, r, vT, rho, k, sigma, lambda, muJ, vJ)
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
  
  cf1 <- om1i * ( r * tau) #ok
  
  cf2 <- vT * k/(sigma^2) * ((k - rho * sigma * om1i - d) * tau - 2 * log((1 - g * exp(-d * tau))/(1 - g)))
  
  cf_jump <- -lambda * muJ * om1i * tau + lambda * tau * ((1 + muJ)^(om1i) * exp(vJ * (om1i/2) * (om1i - 1)) - 1)
  
  A = cf1+cf2+cf_jump
  
  B = 1/sigma^2 * (k - rho * sigma * om1i - d) * (1 - exp(-d * tau))/(1 - g * exp(-d * tau))
  
  w = 2*k/sigma^2
  nu = 2*k*vT/sigma^2
  
  res = exp(A)*(w/(w-B))^nu
  return(res)

}




# adapted from https://stackoverflow.com/questions/10029956/calculating-a-density-from-the-characteristic-function-using-fft-in-r/10038141
characteristic_function_to_density <- function(
  phi, # characteristic function; should be vectorized
  n,   # Number of points, ideally a power of 2
  a, b # Evaluate the density on [a,b[
) {
  i <- 0:(n-1)            # Indices
  dx <- (b-a)/n           # Step size, for the density
  x <- a + i * dx         # Grid, for the density
  dt <- 2*pi / ( n * dx ) # Step size, frequency space
  c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
  d <-  n/2 * dt          # (center the interval on zero)
  t <- c + i * dt         # Grid, frequency space
  phi_t <- phi(t)
  # print(dt)
  X <- exp( -(0+1i) * i * dt * a ) * phi_t
  Y <- fft(X)
  density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
  data.frame(
    i = i,
    t = t,
    characteristic_function = phi_t,
    x = x,
    density = Re(density)
  )
}