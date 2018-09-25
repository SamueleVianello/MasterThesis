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


############## Significativity of correlation ##############



