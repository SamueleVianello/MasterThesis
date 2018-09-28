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
  
  yy= seq(from = 0, to = 0.1,length.out = N+1)
  xx = sqrt( (a*yy^2 - 2*b*yy + c)/d)
  
  if (!full){
    min_sigma = min(xx)
    idx = which(xx>=min_sigma)
    xx= xx[idx]
    yy= yy[idx]
  }

  if(plot){
    plot(xx,yy,type ='l')
    points(diag(sd),r,col='blue')
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
  
  mu = drop((a*expected_return - b)/d)
  lambda = drop( (d - a*b*expected_return + b^2)/(a*d))
  
  
  if (!is.na(expected_return)){
    
    w_opt = invS %*% (e*lambda + mu*r)
    
    sigma_opt = sqrt(t(w_opt)%*%S%*%w_opt)
    sigma = sqrt( (a*expected_return^2 - 2*b*expected_return + c)/d)
    
    # test on the results:
    if (abs(sigma - sigma_opt) > 1e-6){
      stop("computations of minimum standard deviation yield different results.")
    }
    
    print(paste("Minimum standard deviation for given expected return is", sigma_opt))
    return(w_opt)
  }
  else if (!is.na(sd)){
    if(sd < sqrt(1/a)){
      stop("Impossible to obtain such a small risk with given assets.")
    }
    
    expected_return_opt = max(c(b) )
  }

  
}