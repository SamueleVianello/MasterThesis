# Controlling the numerical inversion of the cf to pdf

my_trap_inversion= function(cf, x, N_nodes,M_upper,...){
  # cf : the characteristic function
  # x : the point at which to compute the pdf
  # N_nodes : number of nodes
  # M_upper : upper extremum for the integration interval [0, M_upper]
  
  d_t = M_upper/N_nodes
  t_j = (1:(N_nodes-1))*d_t  # first and last are calculated separately, see next computation
  sum_cos = sum(cos(t_j*x)*Re(cf(t_j,...)))
  sum_sin = sum(sin(t_j*x)*Im(cf(t_j,...)))
  res= 0.5 + 0.5 * cos(M_upper*x)*Re(cf(M_upper,...)) + sum_cos+ sum_sin
  # print(paste("sum_cos=", sum_cos, "sum_sin=", sum_sin))
  return(res*d_t/pi)
}


fourier_cosine = function(cf, x, N,a,b, ...){
  delta = b-a
  k = 0:(N-1)
  F_k = rep(0,N)
  F_k = Re( cf(k*pi/delta, ...)* exp(-1i*k*pi*a/delta))*2/delta
  cos_k = cos(k*pi*(x-a)/delta)
  cos_k[1]=cos_k[1]*0.5
  
  res = sum(F_k*cos_k)
  return(res)
}



source("BatesModel.R")

#### On parameters calibrated on a gaussian
# params = c(0.1534220, 0.9592630, 1.4682972, 0.0000100, 0.8181014)
params = c(-0.0519, 0.488,0.282,0.524,-0.435)
r=params[1]
k= params[2]
eta=params[3]
theta=params[4]
rho=params[5]

# feller condition
print(paste("Feller condition: ", 2*k*eta > theta^2))

x_0 =.32
sigma_0 = 0.5


xx = seq(-20+x_0,20+x_0, length.out=100)
dt = rep(10,length(xx))

yy_trap = xx*0
yy_cos = xx*0
yy_heston = pdfHeston(x=xx,x_0 = x_0, sigma_0 = sigma_0, dt = dt,
                                r = r, k = k, eta = eta, theta = theta, rho = rho)
for(i in 1:length(xx)){
  yy_trap[i]= my_trap_inversion(cf = my_cfHeston, x = xx[i], N_nodes = 100, M_upper = 10, 
                                x_0=x_0, tau=dt[i], r=r, v0=sigma_0^2, vT=eta, rho=rho, k=k, sigma=theta)
  yy_cos[i]= fourier_cosine(cf = my_cfHeston, x = xx[i], a=-10,b=10, N=100,
                               x_0=x_0, tau=dt[i], r=r, v0=sigma_0^2, vT=eta, rho=rho, k=k, sigma=theta)
}

plot(xx,yy_trap)
grid()
lines(xx,yy_heston, col='green')
lines(xx, yy_cos, col='blue')


t1=Sys.time()
yy_heston = pdfHeston(x=xx,x_0 = x_0, sigma_0 = sigma_0, dt = dt,
                      r = r, k = k, eta = eta, theta = theta, rho = rho)
t2=Sys.time()
t2-t1

t1=Sys.time()
for(i in 1:length(xx)){
  yy_trap[i]= my_trap_inversion(cf = my_cfHeston, x = xx[i], N_nodes = 100, M_upper = 10, 
                                x_0=x_0, tau=dt[i], r=r, v0=sigma_0^2, vT=eta, rho=rho, k=k, sigma=theta)
}
t2=Sys.time()
t2-t1

t1=Sys.time()
for(i in 1:length(xx)){
  yy_cos[i]= fourier_cosine(cf = my_cfHeston, x = xx[i], a=-10,b=10, N=100,
                            x_0=x_0, tau=dt[i], r=r, v0=sigma_0^2, vT=eta, rho=rho, k=k, sigma=theta)
}
t2=Sys.time()
t2-t1



###

params = c(-0.734311, 1.00000e-05, 0.0148842,  1.00000,-.99900000)
r=params[1]
k= params[3]
eta=params[2]/params[3]
theta=params[4]
rho=params[5]

# feller condition
print(paste("Feller condition: ", 2*k*eta > theta^2))

cum_returns_eurostoxx = cumsum(rev(eurostoxx))
time_intervals = rev(as.double((as.Date((btc_date)) - as.Date(btc_date[length(btc_date)]) + 1)/365 ))

N = length(cum_returns_eurostoxx)
dn=100

x_0 = log(my_data$EUROSTOXX50[N])
sigma_0 = sd( cum_returns_eurostoxx[1:25])

xx = seq(-20,20, length.out=1000)
dt = rep(time_intervals[1],length(xx))

yy_trap = xx*0
yy_heston = pdfHeston(x=xx,x_0 = x_0, sigma_0 = sigma_0, dt = dt,
                      r = r, k = k, eta = eta, theta = theta, rho = rho)
for(i in 1:length(xx)){
  yy_trap[i]= my_trap_inversion(cf = my_cfHeston, x = xx[i], N_nodes = 1000, M_upper = 10, 
                                x_0=x_0, tau=dt[i], r=r, v0=sigma_0^2, vT=eta, rho=rho, k=k, sigma=theta)
}

plot(xx,yy_trap)
lines(xx,yy_heston, col='green')


