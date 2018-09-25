load("returns.Rda")
load("data.Rda")

library(Hmisc)

attach(my_returns)


############ Plots of coupled returns ####################

x11()
par(mfrow= c(3,5))
plot(btc,bric, pch=20)
plot(btc,sp500, pch=20)
plot(btc,eurostoxx, pch=20)
plot(btc,gold, pch=20)
plot(btc,wti, pch=20)
plot(btc,grain, pch=20)
plot(btc,metal, pch=20)
plot(btc,eur, pch=20)
plot(btc,gbp, pch=20)
plot(btc,chf, pch=20)
plot(btc,jpy, pch=20)
plot(btc,pan_euro, pch=20)
plot(btc,pan_us, pch=20)


# No apparent correlation from graphical inspection


############## Significance of correlation ##############


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

cor(btc,sp500)

PermutationTestCorr(btc,sp500)
rcorr(btc,sp500, type = "pearson")$P['x','y']

# Notice that Spearman is a non parametric test

p_values= data.frame(type = c("Pearson", "Permutation", "Spearman"))
for(i in 2:(dim(my_returns)[2]%/%2)){
  cors = c(PermutationTestCorr(btc,my_returns[,2*i]),
           rcorr(btc,my_returns[,2*i], type = "pearson")$P['x','y'],
           rcorr(btc,my_returns[,2*i], type = "spearman")$P['x','y'])
  p_values[colnames(my_returns)[2*i]] = cors
}

p_values




################## Rolling Correlation ###################################

earliest = min(btc_date)
latest = max(btc_date)














