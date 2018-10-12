# test alternative calibration from 
# Dragulescu and Yakovenko, 


integrand_dragonescu = function(u,x,k,eta,theta, rho, t){
  a1 = 1i*u*x
  G = k+ 1i*rho*theta*u
  omega = sqrt(G^2 + theta^2 * (u^2 -1i*u))
  a2 = k * eta* G* t/theta^2
  a3=-2*k*eta/theta^2 * log( cosh(omega*t*0.5)+ (omega^2 - G^2 + 2*k*G )/(2*k*omega)*sinh(omega*t*0.5))
  
  return(Re(exp(a1+a2+a3)))                          
}



model_pdf = function(x, k,eta,theta, rho, t){
  
  y = rep(0,length(x))
  for(i in 1:length(x)){
    y[i] = integrate(f = integrand_dragonescu, x = x[i], k = k, eta =eta, theta=theta, rho = rho, t = 1/255,
                     lower = -20, upper=20, abs.tol = 1e-5)$value/(2*pi)
  }
  
  return(y)
}


obj = function(params, x_center, empiric_pdf, t){
  
  model_den = model_pdf(x_center, k=params[1], eta=params[2], theta = params[3], rho = params[4],t=t)
  res = sum((model_den-empiric_pdf)^2)
  # print(res)
  # plot(x_center, empiric_pdf)
  # points(x_center, model_den)
  return(res)
}



#########################
# simple plot of integrand
params = c(-0.0519, 0.488,0.282,0.524,-0.435)
r=params[1]
k= params[2]
eta=params[3]
theta=params[4]
rho=params[5]


u = seq(-30,30,length.out = 100)
yy = integrand_dragonescu(u = u,x = 0.6,k=k,eta = eta, theta = theta, rho = rho, t = 1/2)
plot(u,yy)

###############################


attach(my_returns)
mu= mean(eurostoxx)
dt = 1/255
adjusted_returns = eurostoxx - mu


x_min = min(adjusted_returns)
x_max = max(adjusted_returns)

N_bin = 100

# binsize
dx= (x_max-x_min)/N_bin

x_sep = x_min+(0:N_bin)*dx
x_mid = x_sep[1:(N_bin)]+dx/2

out =hist(adjusted_returns, breaks=x_sep, freq = FALSE)

p_mid = out$density




# y[i] = integrate(f = integrand_dragonescu, x = x_centers[i], t = 1/255, eta, rho = rho,k = k,sigma = theta,
#                  lower = lower, upper=upper, abs.tol = 1e-5)$value/pi


# print("Starting calibration using DEoptim...")
# control_list_deoptim = list(itermax = 30, NP = 200, strategy = 1,trace=1)
# 
# start_time_deoptim <- Sys.time()
# outDE <- DEoptim(obj,
#                  lower = c(1e-5,1e-5,1e-5,-1), upper = c(5,1,1,1), control = control_list_deoptim,
#                  x_center=x_mid, t = 1/255, empiric_pdf = p_mid)
# end_time_deoptim <- Sys.time()
# 
# initial=outDE$optim$bestmem



print("Starting calibration using nlminb...")

start_time_nlminb <- Sys.time()

initial = c(runif(4,), runif(1,min=-1))
            
out_nlminb = nlminb(start = initial, objective = obj, lower = c(1e-5, 1e-5, 1e-5,-1), upper = c(5,1,1,1),
                    x_center=x_mid, t = 1/255, empiric_pdf = p_mid,
                    control=list(eval.max = 1000,iter.max = 100, trace = 1))
print(out_nlminb)
end_time_nlminb <- Sys.time()

plot(x_mid, p_mid)
yy= model_pdf(x_mid, k = out_nlminb$par[1],eta = out_nlminb$par[2],theta= out_nlminb$par[3], rho = out_nlminb$par[4], t=1/255)
lines(x_mid,yy)
