plot_simulations = function(times, sims, new_window = TRUE){
  # plots simulations of a stochastic process
  if (new_window) {x11()}
  
  max_s = max(sims)
  min_s = min(sims)
  
  plot(t, sims[1,], type = 'l',col= rainbow(10)[1], ylim=c(min_s - (max_s-min_s)*.1, max_s + (max_s-min_s)*.1))
  for (i in 2:dim(sims)[1]) {
    lines(t, sims[i,], type = 'l',col= rainbow(10)[i%%10 +1])
  }
}