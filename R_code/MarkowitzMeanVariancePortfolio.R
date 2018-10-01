###############################################################
############## MARKOVITZ MEAN VARIANCE PORTFOLIO ##############
###############################################################

EfficientFrontier = function(r,S,full=TRUE,plot=TRUE, N=100){
  # INPUT 
  # r: expected returns of the assets
  # S: covariance matrix of the asset returns
  # full: computes only upper section of frontier if FALSE
  # plot: if TRUE the frontier is plotted
  # N: how many expected returns to take into consideration to plot frontier
  
  
  invS = solve(S)
  e = matrix(rep(1,length(r)),nrow = length(r),ncol = 1) # unit vector 
  
  a = drop(t(e)%*% invS %*% e)
  b = drop(t(e)%*% invS %*% r)
  c = drop(t(r)%*% invS %*% r)
  d = drop(a*c - b^2)
  
  yy= seq(from = 0, to = 0.01,length.out = N+1)
  xx = sqrt( (a*yy^2 - 2*b*yy + c)/d)
  
  if (!full){
    min_sigma = min(xx)
    idx = which(yy>=yy[which(xx==min_sigma)])
    # print(idx)
    xx= xx[idx]
    yy= yy[idx]
  }

  if(plot){
    min_y = min(c(yy,r))
    max_y = max(c(yy,r))
    plot(xx,yy,type ='l', ylim = c(min_y,max_y))
    points(sqrt(diag(S)),r,col='blue',pch='+')
    #legend(c("Efficient Frontier","Single Assets"))
  }
  res = list(sigma = xx, expected_return=yy)
  return(res)
}





OptimalAllocation= function(r,S, expected_return = NA, sd = NA){
  
  if(is.na(expected_return) & is.na(sd)){
    stop("Please specify expected return or standard deviation.")
  }
  
  invS = solve(S)
  e = matrix(rep(1,length(r)),nrow = length(r),ncol = 1) # unit vector 
  
  a = drop(t(e)%*% invS %*% e)
  b = drop(t(e)%*% invS %*% r)
  c = drop(t(r)%*% invS %*% r)
  d = drop(a*c - b^2)
  
  
  if (!is.na(expected_return)){
    mu = drop((a*expected_return - b)/d)
    lambda = drop( (d - a*b*expected_return + b^2)/(a*d))
    
    w_opt = invS %*% (e*lambda + mu*r)
    
    sigma_opt = sqrt(t(w_opt)%*%S%*%w_opt)
    sigma = sqrt( (a*expected_return^2 - 2*b*expected_return + c)/d)
    
    # test on the results:
    if (abs(sigma - sigma_opt) > 1e-6){
      stop("Computations of minimum standard deviation yield different results.")
    }
    
    print(paste("Minimum standard deviation for given expected return is", sigma_opt))

  }
  else if (!is.na(sd)){
    if(sd < sqrt(1/a)){
      stop("Impossible to obtain such a small risk with given assets.")
    }

    expected_return_opt =  (b+sqrt(b^2-a*(c-d*sd^2)))/a
    mu = drop((a*expected_return_opt - b)/d)
    lambda = drop( (d - a*b*expected_return_opt + b^2)/(a*d))
    
    w_opt = invS %*% (e*lambda + mu*r)
    
    sigma_opt = sqrt(t(w_opt)%*%S%*%w_opt)
   
    
    # test on the results:
    if (abs(sd - sigma_opt) > 1e-6){
      stop("Computation of minimum standard deviation yields a different result.")
    }
    print(paste("Maximum expected return for given standard deviation is",expected_return_opt))
  }
  
  return(w_opt)
}