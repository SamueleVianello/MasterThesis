library(NMOF)


u = 1

r=0.8
t=1
S_0 = 1
sigma_0 = 0.05
theta =0.8
rho=0.3
k=2
eta = 0.5
lambda=0.5
mu_j = -0.1
sigma_j = 0.17


a = 1i*rho*theta*u

g = sqrt(theta^2 * (u^2+1i*u)+ (k-a)^2)

# phi_j = exp( t*lambda*(exp(-sigma_j^2*u^2*0.5+1i*u*(log(1+mu_j)-sigma_j^2*0.5))-1))
# 
# phi_d = exp(k*eta*t*(k-a)/theta^2 + 1i*u*t*(r-lambda*mu_j)+1i*u*log(S_0))  *  exp(-(u^2+1i*u)*sigma_0/(g*1/tanh(g*t*0.5)+k-a))  / (cosh(g*t*0.5)+(k-a)*sinh(g*t*0.5)/g)^(2*k*eta/theta^2)
# 
# 
# 
# phi_d*phi_j
# cfBates(om = u,S = S_0,tau = t,r = r,q = 0,v0 = sigma_0,vT = eta,rho = rho,k = k,
#         sigma = theta,lambda = lambda,muJ = mu_j,vJ = sigma_j)

e1 = exp()

h = (k-a-g)/(k-a+g)
ex2= exp(k*eta/theta^2 * ((k-a -g)*t-2*log((1-h*exp(-g*t))/(1-h))))
ex3= exp(sigma_0^2/theta^2*(k-a-g)*(1-exp(-g*t))/(1-h*exp(g*t)))

e1*ex2*ex3

cfHeston(om = u,S = S_0,tau = t,r = r, q = 0,v0 = sigma_0^2,vT = eta,rho = rho,k = k,sigma = eta)
         #         sigma = theta,lambda = lambda,muJ = mu_j,vJ = sigma_j)



f=function(x, S,tau ,r , q ,v0,vT,rho,k,sigma){
  return(Re(cfHeston(om = x,S = S_0,tau = t,r = r, q = 0,v0 = sigma_0^2,vT = eta,rho = rho,k = k,sigma = eta)
            * exp(-1i*log(S)*x)))
}

integrate(f = f, S = 1,tau = t,r = r, q = 0,v0 = sigma_0^2,vT = eta,rho = rho,k = k,sigma = eta,
          lower = 0, upper=50)







my_cfHeston=function (om, z, tau, r, v0, vT, rho, k, sigma)
{
  if (sigma < 1e-08) 
    sigma <- 1e-08
  d <- sqrt((rho * sigma * (0+1i) * om - k)^2 + sigma^2 * ((0+1i) * om + om^2))
  
  g <- (k - rho * sigma * (0+1i) * om - d)/(k - rho * sigma * (0+1i) * om + d)
  
  cf1 <- (0+1i) * om * (z + (r) * tau) 
  
  cf2 <- vT * k/(sigma^2) * ((k - rho * sigma * (0+1i) * om - d) * tau - 2 * log((1 - g * exp(-d * tau))/(1 - g)))
  
  cf3 <- v0/sigma^2 * (k - rho * sigma * (0+1i) * om - d) * (1 - exp(-d * tau))/(1 - g * exp(-d * tau))
  exp(cf1 + cf2 + cf3)
}

f2=function(x, z,tau ,r , q ,v0,vT,rho,k,sigma){
  return(Re(my_cfHeston(om = x,z = z,tau = t,r = r,v0 = sigma_0^2,vT = eta,rho = rho,k = k,sigma = eta)
            * exp(-1i*z*x)))
}


xx = seq(from=-10, to = 10, by =20/1000)
yy=rep(0,length(xx))
yy2=rep(0,length(xx))
for( i in 1:length(xx)){
  yy[i] = integrate(f = f, S = exp(xx[i]),tau = t,r = r, q = 0,v0 = sigma_0^2,vT = eta,rho = rho,k = k,sigma = eta,
                    lower = 0, upper=50)$value/pi
  # yy2[i]= integrate(f = f2, z = (xx[i]),tau = t,r = r,v0 = sigma_0^2,vT = eta,rho = rho,k = k,sigma = eta,
  #                    lower = 0, upper=50)$value/pi
}
plot(xx,yy,type='l')
lines(xx,yy2,col='blue')
sum(yy*(xx[2]-xx[1]))
