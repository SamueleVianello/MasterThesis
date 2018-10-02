load("returns.Rda")
load("data.Rda")
load("results.Rda")

source("MultivariateMertonModel.R")
source("CalibrationMVMerton.R")

averaged_mu = colSums(results$full_mu)/6
averaged_delta = colSums(results$full_delta)/6
averaged_theta = colSums(results$full_theta)/6
averaged_lambda = colSums(results$full_lambda)/6


attach(my_returns)
calib_nojump = CalibrateMVMerton(x=btc, n=1, dt = 1/255, min_theta = 0, max_theta = 0, min_delta = 0, max_delta = 0)

calib_pos = CalibrateMVMerton(x=btc, n=1, dt = 1/255, min_theta =0.1, max_theta = 1 )
calib_neg = CalibrateMVMerton(x=btc, n=1, dt = 1/255, min_theta = -1, max_theta = -0.1 )


calib_unconstr =  CalibrateMVMerton(x=btc, n=1, dt = 1/255)



xx = seq(from = -5, to = 5, by = 0.001)
# yy_averaged = MultivariateMertonPdf_1asset_nocommon(x = xx,dt = 1/255,mu = averaged_mu[1],S = results$covariance[1,1],
#                                            theta = averaged_theta[1],delta = averaged_delta[1],lambda = averaged_lambda[1])
yy_pos = MultivariateMertonPdf_1asset_nocommon(x = xx, dt = 1/255, mu = calib_pos$mu, S = calib_pos$S,
                                           theta = calib_pos$theta, delta = calib_pos$delta, lambda = calib_pos$lambda)
yy_neg = MultivariateMertonPdf_1asset_nocommon(x = xx, dt = 1/255, mu = calib_neg$mu, S = calib_neg$S,
                                               theta = calib_neg$theta, delta = calib_neg$delta, lambda = calib_neg$lambda)
yy_unconstr = MultivariateMertonPdf_1asset_nocommon(x = xx, dt = 1/255, mu = calib_unconstr$mu, S = calib_unconstr$S,
                                               theta = calib_unconstr$theta, delta = calib_unconstr$delta, lambda = calib_unconstr$lambda)


h =hist(x = btc, breaks = 50, freq = FALSE)
curve(dnorm(x, mean=mean(btc), sd=sd(btc)), add=TRUE,col="black")
lines(xx, yy_pos,col='green')
lines(xx, yy_neg,col='red')
lines(xx, yy_unconstr,col='blue')

legend("topleft", legend = c("Black-Scholes","Merton_pos","Merton_neg", "Merton_unconstr"),col=c('black','green','red', 'blue'),
       lwd=3,lty=c(1,1,1),cex=0.75)

likelihood = rbind(calib_neg$objective_function,calib_pos$objective_function,calib_unconstr$objective_function)
rownames(likelihood)=c("neg", "pos", "unconstr")
likelihood

lambdas = rbind(calib_neg$lambda,calib_pos$lambda,calib_unconstr$lambda)
rownames(lambdas)=c("neg", "pos", "unconstr")
lambdas


