# empirical VaR portfolio optimization

source("MarkowitzMeanVariancePortfolio.R")
load("returns.Rda")
load("data.Rda")
load("results.Rda")

library(nlminb)
library(nloptr)
library(Rsolnp)
library(GenSA)


N_samples=dim(my_returns)[1]
asset_names = colnames(my_returns[,2*(1:14)])


# rendimenti percentuali:
percentage_returns = exp(my_returns[,2*(1:14)])

 
hist(percentage_returns$sp500, freq = FALSE)


w = rep(1/14,14)

mean(ptf_daily)^255

# annual percentage simulated return vector

N_assets = length(w)

N_sim= 1000 # takes almost 5 mins for 1000!

sim_mat = matrix(rep(0,N_sim*N_assets), nrow = N_sim)
x=rep(1,N_assets)

t_beg = Sys.time()
for(k in 1:N_sim){
  x = x=rep(1,N_assets)
  samples = sample(N_samples,255, replace = TRUE)
  for( i in 1:255){
    x = x* percentage_returns[samples[i],]
  }
  sim_mat[k,]=as.matrix(x)
}
t_end = Sys.time()

t_end-t_beg

annual_returns = sim_mat %*% w

hist(log(annual_returns), breaks = 80, main = "Naif allocation annual log-return")
abline(v = quantile(log(annual_returns),probs = 0.05), col = 'blue')

VaR = quantile(annual_returns,probs = c(0.01,0.05,0.1))
print(paste("percentage VaR for naif allocation at",  c(0.01,0.05,0.1)*100 , "% level is",VaR))


ptf_var = function(w,  sims, alpha, target_return){
  annual = sims%*%w
  res = quantile(x = annual, probs = alpha)
  return(res)
}



# constr_matrix = rep(1,N_assets)
# constr_matrix = rbind(constr_matrix, rep(-1,N_assets))
# 
# mean_asset_return = t(sim_mat)%*% rep(1, N_sim)/N_sim
# constr_matrix = rbind(constr_matrix, t(mean_asset_return))
# constr_matrix = rbind(constr_matrix, diag(rep(1,N_assets)))
# 
# target = 1.10
# b = c(1,-1,target, rep(0,N_assets))-1e-8
# 

# initial = rep(1/N_assets,N_assets)


# # Not working
# sol = constrOptim(theta = initial,f = ptf_var,
#                   ui = constr_matrix,ci = b, 
#                   sims = sim_mat, alpha = 0.05, 
#                   grad = NULL, method = "Nelder-Mead", control = list(trace=1))
# sol

cbind(constr_matrix %*% initial,b)



f_ineq = function(w, sims, alpha,target_return){
  c1 = sum(w)-1
  c2 = sum(-w)+1
  N_s = dim(sims)[1]
  c3 = sum(t(sims)%*% rep(1, N_s)*w) /N_s - target_return
  
  return(c(c1, c2, c3))
}

# initial =c(0.5,0.5, rep(0,N_assets-2))
initial = rep(1/N_assets,N_assets)
# initial[1] =2/14
# initial[2] = 0
# Also not working: doesn't move from initial point...
res = cobyla(x0=initial, fn = ptf_var, lower = rep(0,N_assets), hin = f_ineq, # nl.info = FALSE,
               sims = sim_mat, alpha = 0.05, target_return = 1.10,
               nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000))

w_opt=matrix(res$par, nrow = N_assets)
rownames(w_opt)=asset_names
w_opt

mean(sim_mat%*%(w_opt))
ptf_var(w = (w_opt), sims =sim_mat, alpha = 0.05,target_return = 1.2)

# everything on oil wti!?

#  *** Global optimization with constraints included in obj function ***

ptf_var_w_constr= function(w, sims, alpha,target_return){
  toll = 1e-7
  N_s = dim(sims)[1]
  if ( sum(t(sims)%*% rep(1, N_s)*w) /N_s - target_return < 0 ){
    return(1)
    }
  else if(sum(w)>1+toll || sum(w)<1-0.001){
    return(1)
  }
  else {
    annual = sims%*%w
    res = quantile(x = annual, probs = alpha)
    return(res)
  }
}

initial =c(1, rep(0,N_assets-1))
# initial = rep(1/N_assets,N_assets)
sol = GenSA(par = initial, fn = ptf_var_w_constr,lower = rep(0,N_assets), upper = rep(1,N_assets),
            sims = sim_mat, alpha = 0.05, target_return = 1.20,
            control = c(verbose=TRUE, temperature = 10000))
sol$par
sum(sol$par)
