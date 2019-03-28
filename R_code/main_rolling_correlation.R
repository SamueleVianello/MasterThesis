library(mondate)
load("returns.Rda")
load("data.Rda")
attach(my_returns)



# Funtcion to perform permutation test ########################################
PermutationTestCorr = function(x,y=0, N=2000){
  # Two sided permutation test for correlation
  # H0: rho = 0   vs    H1: rho!=0
  
  if (!is.null(dim(x)[2])){
    if ( dim(x)[2] == 2){
      y=x[,2]
      x=x[,1]
    }
  }
  else if (length(x)!=length(y)){
    stop("Error: samples should have same size.")
  }
  
  n = length(x)
  
  r_sample= cor(x,y)
  
  larger = 0
  for(i in 1:N){
    y_perm = y[sample(n,n)]
    r_perm = cor(x,y_perm)
    if (abs(r_sample)<abs(r_perm))
      larger = larger+1
  }
  
  # p-value is the percentage of r_perm absolutely greater than r_sample
  p = larger/N
  return(p)
} 
############################################################################


dates = as.Date(btc_date, origin="1899-12-30")
earliest = min(dates)
latest = max(dates)

w = 36 # rolling window in months
step = 1 # moving step in months

d = as.Date(mondate(latest, timeunits = )-step,origin="1899-12-30")
beg_times = as.Date(mondate(latest)-w)
end_times = latest
while( as.Date(mondate(d)-w) > earliest){
  beg_times = c(as.Date(mondate(d)-w), beg_times)
  end_times = c( d, end_times)
  d = as.Date(mondate(d)-step)
}

rbind(beg_times,end_times)

n_wind = length(beg_times)
n_assets = (dim(my_returns)[2]%/%2)


roll_corr = matrix(rep(NA,n_wind*(n_assets-1)),n_wind,n_assets-1)
roll_pvalue = matrix(rep(NA,n_wind*(n_assets-1)),n_wind,n_assets-1)

for (i in 1:n_wind){
  idx = which(dates >= beg_times[i] & dates<=end_times[i])
  print(idx[1:5])
  correlations = cor(btc[idx],my_returns[idx,2*(2:n_assets)])
  roll_corr[i,]=correlations
  
  p_values= rep(NA,n_assets-1)
  for(j in 2:n_assets){
    # # Using Pearsons
    # p_values[j-1] = rcorr(btc[idx],my_returns[idx,2*j], type = "pearson")$P['x','y']
    p_values[j-1] = PermutationTestCorr(btc[idx],my_returns[idx,2*j])
  }
  roll_pvalue[i,] = p_values
}

colnames(roll_corr)= colnames(my_returns[2*(2:(n_assets))])
colnames(roll_pvalue)= colnames(my_returns[2*(2:(n_assets))])

roll_corr_3y = roll_corr
roll_pvalue_3y = roll_pvalue
end_times_3y = end_times
beg_times_3y = beg_times




#########

w = 18 # rolling window in months
step = 1 # moving step in months

d = as.Date(mondate(latest)-step)
beg_times = as.Date(mondate(latest)-w)
end_times = latest
while( as.Date(mondate(d)-w) > earliest){
  beg_times = c(as.Date(mondate(d)-w), beg_times)
  end_times = c( d, end_times)
  d = as.Date(mondate(d)-step)
}

# rbind(beg_times,end_times)

n_wind = length(beg_times)
n_assets = (dim(my_returns)[2]%/%2)


roll_corr = matrix(rep(NA,n_wind*(n_assets-1)),n_wind,n_assets-1)
roll_pvalue = matrix(rep(NA,n_wind*(n_assets-1)),n_wind,n_assets-1)

for (i in 1:n_wind){
  idx = which(dates >= beg_times[i] & dates<=end_times[i])
  print(idx[1:5])
  correlations = cor(btc[idx],my_returns[idx,2*(2:n_assets)])
  roll_corr[i,]=correlations
  
  p_values= rep(NA,n_assets-1)
  for(j in 2:n_assets){
    # # Using Pearsons
    # p_values[j-1] = rcorr(btc[idx],my_returns[idx,2*j], type = "pearson")$P['x','y']
    p_values[j-1] = PermutationTestCorr(btc[idx],my_returns[idx,2*j])
  }
  roll_pvalue[i,] = p_values
}


colnames(roll_corr)= colnames(my_returns[2*(2:(n_assets))])
colnames(roll_pvalue)= colnames(my_returns[2*(2:(n_assets))])

# p_bonf_18m = roll_pvalue*0
# for (i in 1:ncol(p_bonf)) {
#   p_bonf_18m[,i]=p.adjust(roll_pvalue[,i],method = "bonf")
# }


library(zoo)
zoo_roll_corr = zoo(roll_corr, order.by = end_times)

########
require(DescTools)

y_min = -0.4
y_max =  0.4

# plot against stock indices
x11()
layout(matrix(c(1,2,3,4,1,2,3,4,5,6,7,8),nrow= 3,ncol=4, byrow=TRUE))
plot(end_times, roll_corr[,'bric'], type = 'l', ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'bric'])
lines(end_times_3y, roll_corr_3y[,'bric'], type = 'l', ylab = "Bric",col='blue')
title("bric")
grid()

plot(end_times, roll_corr[,'sp500'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'sp500'])
lines(end_times_3y, roll_corr_3y[,'sp500'], type = 'l', ylab = "sp500",col='blue')
title("sp500")
grid()

plot(end_times, roll_corr[,'eurostoxx'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'eurostoxx'])
lines(end_times_3y, roll_corr_3y[,'eurostoxx'], type = 'l', ylab = "eurostoxx",col='blue')
title("eurostoxx")
grid()

plot(end_times, roll_corr[,'nasdaq'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'eurostoxx'])
lines(end_times_3y, roll_corr_3y[,'nasdaq'], type = 'l', ylab = "nasdaq",col='blue')
title("nasdaq")
grid()


Mar(top=1)
plot(end_times, roll_pvalue[,'bric'],type = 'b',pch=16,ylim = c(0,1), ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'bric'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')


plot(end_times, roll_pvalue[,'sp500'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'sp500'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')

plot(end_times, roll_pvalue[,'eurostoxx'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'eurostoxx'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')

plot(end_times, roll_pvalue[,'nasdaq'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'nasdaq'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')

# dev.copy2pdf(file="rolling_stocks.pdf")
# dev.off()



# plot against commodities
x11()
layout(matrix(c(1,1,5,2,2,6,3,3,7,4,4,8),nrow= 3,ncol=4, byrow=FALSE))
plot(end_times, roll_corr[,'gold'], type = 'l', ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'bric'])
lines(end_times_3y, roll_corr_3y[,'gold'], type = 'l', ylab = "gold",col='blue')
title("gold")
abline(v=axis.Date(1, x=pretty(index(zoo_roll_corr))),col = "lightgray", lty = "dotted", lwd = par("lwd"))
grid(nx=NA,  ny =NULL)

plot(end_times, roll_corr[,'wti'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'sp500'])
lines(end_times_3y, roll_corr_3y[,'wti'], type = 'l', ylab = "wti",col='blue')
title("wti")
abline(v=axis.Date(1, x=pretty(index(zoo_roll_corr))),col = "lightgray", lty = "dotted", lwd = par("lwd"))
grid(nx=NA,  ny =NULL)

plot(end_times, roll_corr[,'grain'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'eurostoxx'])
lines(end_times_3y, roll_corr_3y[,'grain'], type = 'l', ylab = "grain",col='blue')
title("grain")
abline(v=axis.Date(1, x=pretty(index(zoo_roll_corr))),col = "lightgray", lty = "dotted", lwd = par("lwd"))
grid(nx=NA,  ny =NULL)

plot(end_times, roll_corr[,'metal'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'eurostoxx'])
lines(end_times_3y, roll_corr_3y[,'metal'], type = 'l', ylab = "metal",col='blue')
title("metal")
abline(v=axis.Date(1, x=pretty(index(zoo_roll_corr))),col = "lightgray", lty = "dotted", lwd = par("lwd"))
grid(nx=NA,  ny =NULL)


Mar(top=1)
plot(end_times, roll_pvalue[,'gold'],type = 'b',pch=16,ylim = c(0,1), ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'gold'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')
abline(v=axis.Date(1, x=pretty(index(zoo_roll_corr))),col = "lightgray", lty = "dotted", lwd = par("lwd"))

plot(end_times, roll_pvalue[,'wti'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'wti'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')
abline(v=axis.Date(1, x=pretty(index(zoo_roll_corr))),col = "lightgray", lty = "dotted", lwd = par("lwd"))

plot(end_times, roll_pvalue[,'grain'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'grain'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')
abline(v=axis.Date(1, x=pretty(index(zoo_roll_corr))),col = "lightgray", lty = "dotted", lwd = par("lwd"))

plot(end_times, roll_pvalue[,'metal'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'metal'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')
abline(v=axis.Date(1, x=pretty(index(zoo_roll_corr))),col = "lightgray", lty = "dotted", lwd = par("lwd"))

# dev.copy2pdf(file="rolling_commodities.pdf")
# dev.off()

# plot against fx
x11()
# layout(matrix(c(1,1,5,2,2,6,3,3,7,4,4,8),nrow= 3,ncol=4, byrow=FALSE))
layout(matrix(c(1,1,5,2,2,6,3,3,7,4,4,8),nrow= 3,ncol=4, byrow=FALSE), heights = rep(1)/3)
#par(mar=c(0, 0, 1, 1) + 0.0)
plot(end_times, roll_corr[,'eur'], type = 'l', ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'bric'])
lines(end_times_3y, roll_corr_3y[,'eur'], type = 'l', ylab = "eur",col='blue')
title("eur")
grid()

plot(end_times, roll_corr[,'gbp'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'sp500'])
lines(end_times_3y, roll_corr_3y[,'gbp'], type = 'l', ylab = "gbp",col='blue')
title("gbp")
grid()

plot(end_times, roll_corr[,'chf'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'eurostoxx'])
lines(end_times_3y, roll_corr_3y[,'chf'], type = 'l', ylab = "chf",col='blue')
title("chf")
grid()

plot(end_times, roll_corr[,'jpy'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'eurostoxx'])
lines(end_times_3y, roll_corr_3y[,'jpy'], type = 'l', ylab = "jpy",col='blue')
title("jpy")
grid()


Mar(top=1)
plot(end_times, roll_pvalue[,'eur'],type = 'b',pch=16,ylim = c(0,1), ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'eur'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')

plot(end_times, roll_pvalue[,'gbp'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'gbp'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')

plot(end_times, roll_pvalue[,'chf'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'chf'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')

plot(end_times, roll_pvalue[,'jpy'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'jpy'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')

# dev.copy2pdf(file="rolling_fx.pdf")
# dev.off()



# plot against bond indices
x11()
layout(matrix(c(1,2,3,1,2,3,4,5,6),nrow= 3,ncol=3, byrow=TRUE))
plot(end_times, roll_corr[,'bond_europe'], type = 'l', ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'bric'])
lines(end_times_3y, roll_corr_3y[,'bond_europe'], type = 'l', ylab = "bond_europe",col='blue')
title("bond_europe")
grid()

plot(end_times, roll_corr[,'bond_us'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'sp500'])
lines(end_times_3y, roll_corr_3y[,'bond_us'], type = 'l', ylab = "bond_us",col='blue')
title("bond_us")
grid()

plot(end_times, roll_corr[,'bond_eur'], type = 'l',ylab = "roll_corr", ylim = c(y_min,y_max))
#points(end_times, roll_corr[,'sp500'])
lines(end_times_3y, roll_corr_3y[,'bond_eur'], type = 'l', ylab = "bond_eur",col='blue')
title("bond_eur")
grid()

# plot(end_times, roll_corr[,'vix'], type = 'l',ylab = "vix", ylim = c(y_min,y_max))
# #points(end_times, roll_corr[,'sp500'])
# lines(end_times_3y, roll_corr_3y[,'vix'], type = 'l', ylab = "vix",col='blue')
# title("vix")
# grid()



Mar(top=0.1)
plot(end_times, roll_pvalue[,'bond_europe'],type = 'b',pch=16,ylim = c(0,1), ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'bond_europe'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')

plot(end_times, roll_pvalue[,'bond_us'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'bond_us'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')

plot(end_times, roll_pvalue[,'bond_eur'],type = 'b',pch=16,ylim = c(0,1),ylab = "pvalue")
lines(end_times_3y,roll_pvalue_3y[,'bond_eur'],type = 'b',pch=16,col = 'blue')
abline(h = .05,col = 'grey')


# plot(end_times, roll_pvalue[,'vix'],type = 'b',pch=16,ylim = c(0,1),ylab = "vix")
# lines(end_times_3y,roll_pvalue_3y[,'vix'],type = 'b',pch=16,col = 'blue')
# abline(h = .05,col = 'grey')

dev.copy2pdf(file="rolling_bonds.pdf")

dev.off()

