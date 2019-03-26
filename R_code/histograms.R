# histogram comparison of merton, heston and bates



result_H
result_B

dt = 1/250

for (i in 1:16){
  returns= my_returns[,2*i]
  x11()
  hist(returns, freq = FALSE, breaks = 50, main = asset_names[i], xlab = "log-returns", ylab = "density")
  x_min = min(returns)*1.2
  x_max = max(returns)*1.2
  xx = seq(length.out = 1000, from = x_min, to=x_max)
  yy_gauss = dnorm(xx, mean = mean(returns), sd = sd(returns))
  
  yy_merton = MultivariateMertonPdf_1asset_nocommon(x=xx, dt = dt, mu = full_results[i,1], S = full_results[i,2]^2,
                                                    theta = full_results[i,3], delta = full_results[i,4], lambda = full_results[i,5])
  yy_heston = pdfHeston_fft(x=xx,dt = dt,r = result_H[i,1], k = result_H[i,2], eta = result_H[i,3],
                            theta = result_H[i,4],rho = result_H[i,5])
  yy_bates = pdfBates_fft(x=xx,dt = dt,r = result_B[i,1], k = result_B[i,2], eta = result_B[i,3],
                          theta = result_B[i,4],rho = result_B[i,5], mu_j = result_B[i,6],sigma_j =  result_B[i,7],lambda = result_B[i,8])
  
  lines(xx,yy_gauss, col="black",lwd=2)
  lines(xx,yy_merton, col="green",lwd=2)
  lines(xx,yy_heston, col="red",lwd=2)
  lines(xx,yy_bates,col="blue",lwd=2)
  
  legend("toprigh", legend = c("gaussian","merton","heston","bates"),
         col=c("black","green","red","blue"), lwd = 2, lty = 1, cex=0.9, bg = 'white')
  
  
  dev.copy2pdf(file = paste0("hist_",asset_names[i],".pdf"), height = 7, width=7 )
}
