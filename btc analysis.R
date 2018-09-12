
library(quantmod)




btc = getSymbols("BTC-USD",from="2015-08-31",to="2018-08-31", auto.assign = FALSE)

btc =  na.omit(btc[,3])









