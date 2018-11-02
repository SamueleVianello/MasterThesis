# test nlopt usage

feller_constr = function(param,x,x_0,dt,sigma_0 , model,  check_feller, unconditional){
  x = param[2]*param[3]*2/param[4]^2 - 1
  return(-x)
}

constr_grad = function(param){
  x = rep(0,length(param))
  x[2] = 2*param[3]/param[4]^2
  x[3] = 2*param[2]/param[4]^2
  x[4] = -4*param[2]*param[3]/param[4]^3
  return(matrix(-x, nrow = 1))
}

bounds = BoundsCreator(model = "heston",sigma_param = F)
options = list(algorithm = "NLOPT_GN_ISRES", print_level = 2, maxeval=500)


asset = sp500
asset_date = eur_date

cum_returns_asset = cumsum(rev(asset))
time_intervals = rev(as.double((as.Date((asset_date)) - as.Date(asset_date[length(asset_date)]) + 1)/365 ))


N = length(cum_returns_asset)
dn=255*2

# x_0 = log(my_data$EUROSTOXX50[N])
x_0 = 0
sigma_0 = sd( rev(asset)[1:dn])*sqrt(255)

# initial_mu = mean(cum_returns_eur[1:dn]/time_intervals[1:dn])
initial_mu = mean(rev(asset)[1:dn])*255
initial =  c(initial_mu, 0.70, sigma_0^2, 0.1, -0.3)
initial[2]*initial[3]*2 > initial[4]^2


plot(time_intervals[1:dn], cum_returns_asset[1:dn], type='l')
set.seed(1234)

t1=Sys.time()

res =  nloptr(x0=initial, eval_f = negloglikHeston, lb = bounds$lower, ub =bounds$upper,eval_g_ineq = feller_constr,eval_jac_g_ineq = ,opts = options,
              x=cum_returns_asset[1:dn],x_0 = x_0, sigma_0=0, dt = time_intervals[1:dn], model="heston",  check_feller = F, unconditional =T)
