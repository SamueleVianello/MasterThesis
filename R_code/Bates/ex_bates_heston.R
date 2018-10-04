my_cfHeston=
function (om, x_0, tau, r, q, v0, vT, rho, k, sigma) 
{
  if (sigma < 1e-08) 
    sigma <- 1e-08
  d <- sqrt((rho * sigma * (0+1i) * om - k)^2 + sigma^2 * ((0+1i) *om + om^2))
  g <- (k - rho * sigma * (0+1i) * om - d)/(k - rho * sigma * (0+1i) * om + d)
  
  cf1 <- (0+1i) * om * (x_0 + r * tau)
  
  cf2 <- vT * k/(sigma^2) * ((k - rho * sigma * (0+1i) * om - d) * tau - 2 * log((1 - g * exp(-d * tau))/(1 - g)))
  
  cf3 <- v0/sigma^2 * (k - rho * sigma * (0+1i) * om - d) * (1 - exp(-d * tau))/(1 - g * exp(-d * tau))
  
  # print(paste("exp1=", cf1, ", exp2=",cf2, ", exp3=",cf3))
  return(exp(cf1 + cf2 + cf3))
}





my_cfBates=
function (om, x_0, tau, r, q, v0, vT, rho, k, sigma, lambda, muJ, 
          vJ) 
{
  if (sigma < 1e-08) 
    sigma <- 1e-08
  sigma <- max(sigma, 1e-04)
  om1i <- om * (0+1i)
  d <- sqrt((rho * sigma * om1i - k)^2 + sigma^2 * (om1i + om^2))
  
  g <- (k - rho * sigma * om1i - d)/(k - rho * sigma * om1i + d)
  
  cf1 <- om1i * (x_0 + r * tau) #ok
  
  cf2 <- vT * k/(sigma^2) * ((k - rho * sigma * om1i - d) * tau - 2 * log((1 - g * exp(-d * tau))/(1 - g)))
  
  cf3 <- v0/sigma^2 * (k - rho * sigma * om1i - d) * (1 - exp(-d * tau))/(1 - g * exp(-d * tau))
  
  cf4 <- -lambda * muJ * om1i * tau + lambda * tau * ((1 + muJ)^(om1i) * exp(vJ * (om1i/2) * (om1i - 1)) - 1)
  
  #print(paste("exp1=", cf1, ", exp2=",cf2, ", exp3=",cf3, ", exp4= ",cf4))
  
  
  
  res= exp(cf1 + cf2 + cf3 + cf4)
  is.na(res)
  return(res)
}