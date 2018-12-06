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
plot(btc,bond_europe, pch=20)
plot(btc,bond_us, pch=20)
plot(btc,bond_eur, pch=20)


# No apparent correlation from graphical inspection

#**************************************** 
#** VISUAL COMPARISON OF CORRELATIONS  **
#****************************************
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)


load("results.Rda")


attach(results)
colnames(correlation) = colnames(my_returns[,2*(1:17)])
rownames(correlation) = colnames(my_returns[,2*(1:17)])

correl_reverse = correlation[dim(correlation)[1]:1,]

longData<-melt(correl_reverse)
longData<-longData[longData$value!=0,]

p1=ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low="blue", high="red",mid = "white", limits = c(-1,1)) +
  labs( title="Model Correlation") +
  theme_bw() + theme(axis.text.x=element_text(size=12, angle=90,hjust = 0),
                     axis.text.y=element_text(size=12),
                     plot.title=element_text(size=15),
                     axis.title = element_blank() )+
  scale_x_discrete(position = "top")


sample_corr = cor(my_returns[, 2*(1:17)])
sample_cor_reverse = sample_corr[dim(sample_corr)[1]:1,]

longData<-melt(sample_cor_reverse)
longData<-longData[longData$value!=0,]

p2=ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low="blue", high="red",mid = "white", limits=c(-1,1)) +
  labs( title="Sample Correlation") +
  theme_bw() + theme(axis.text.x=element_text(size=12, angle=90,hjust = 0),
                     axis.text.y=element_text(size=12),
                     plot.title=element_text(size=15),
                     axis.title = element_blank() )+
  scale_x_discrete(position = "top")

x11()
grid.arrange(p1,p2,nrow=1)



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


p_values= data.frame(type = c("Correlation","Pearson", "Permutation", "Spearman"))
for(i in 2:(dim(my_returns)[2]%/%2)){
  cors = c(cor(btc,my_returns[,2*i]),
           PermutationTestCorr(btc,my_returns[,2*i]),
           rcorr(btc,my_returns[,2*i], type = "pearson")$P['x','y'],
           rcorr(btc,my_returns[,2*i], type = "spearman")$P['x','y'])
  p_values[colnames(my_returns)[2*i]] = cors
}

p_values















