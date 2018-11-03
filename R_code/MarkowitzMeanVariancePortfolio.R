###############################################################
############## MARKOVITZ MEAN VARIANCE PORTFOLIO ##############
###############################################################

# Simple wrapper for frontier
EfficientFrontier=function(r,S, no_short_sales=NA, N=100,full=FALSE,plot=FALSE,min_r=NA, max_r=NA){
  
  if (is.na(no_short_sales[1])){
    return(EfficientFrontier_unconstr(r=r,S=S,full=full,plot=plot, N=N, max_r = max_r, min_r=min_r))

  }
  else{
    return(EfficientFrontier_constr(r=r,S=S,full=full,plot=plot, N=N, no_short_sales = no_short_sales,
                                    min_r=min_r, max_r = max_r))
  }
  
}

#Simple wrapper for optimal allocation
OptimalAllocation=function(r,S, target_return = NA, sd = NA, no_short_sales=NA){
  if (is.na(no_short_sales[1])){
    return(OptimalAllocation_unconstr(r=r,S=S, target_return = target_return, sd=sd))
  }
  else{
    return(OptimalAllocation_constr(r=r,S=S, target_return = target_return,no_short_sales = no_short_sales))
    
  }
}



OptimalAllocation_unconstr= function(r,S, target_return = NA, sd = NA){
  
  if(is.na(target_return) & is.na(sd)){
    stop("Please specify expected return or standard deviation.")
  }
  
  invS = solve(S)
  e = matrix(rep(1,length(r)),nrow = length(r),ncol = 1) # unit vector
  
  a = drop(t(e)%*% invS %*% e)
  b = drop(t(e)%*% invS %*% r)
  c = drop(t(r)%*% invS %*% r)
  d = drop(a*c - b^2)
  
  
  if (!is.na(target_return)){
    mu = drop((a*target_return - b)/d)
    lambda = drop( (d - a*b*target_return + b^2)/(a*d))
    
    w_opt = invS %*% (e*lambda + mu*r)
    
    sigma_opt = sqrt(t(w_opt)%*%S%*%w_opt)
    sigma = sqrt( (a*target_return^2 - 2*b*target_return + c)/d)
    
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



OptimalAllocation_constr= function(r,S, target_return = NA, no_short_sales){
  
  if(is.na(target_return)){
    stop("Please specify expected return.")
  }

  if(target_return > max(r)){
    warning("Cannot produce such a high return without short-selling.")
    target_return=max(r)-1e-5
  }
  
  library(quadprog)
  
  # number of assets
  n= length(r)
  
  D = 2*S       #times 2 because there is a 1/2 in the implicit formulation
  d = rep(0,n)  # zeros
  
  
  # Constraint on returns
  A = t(r)
  
  b = target_return
  
  # Constraint on weights: sum(w)=1
  A = rbind(A,rep(1,n))
  b= rbind(b,1)
  
  
  # Constraint on short selling ony on given assets
  for (i in no_short_sales){
    A_i = rep(0,n)
    A_i[i]=1
    A = rbind(A,A_i)
    b = rbind(b,0)
  }
  
  sol = solve.QP(Dmat = D, dvec = (d), Amat = t(A), bvec = t(b), meq = 2)
  w_opt= sol$solution
  rownames(w_opt)=rownames(r)
  return(w_opt)
}





EfficientFrontier_unconstr = function(r,S,full=FALSE,plot=FALSE, N=100, max_r=NA, min_r=NA){
  # INPUT
  # r: expected returns of the assets
  # S: covariance matrix of the asset returns
  # full: computes only upper section of frontier if FALSE
  # plot: if TRUE the frontier is plotted
  # N: how many expected returns to take into consideration to plot frontier
  if (is.na(max_r)){
    max_r = max(r)
  }
  if (is.na(min_r)){
    min_r = min(r)
  }

  invS = solve(S)
  e = matrix(rep(1,length(r)),nrow = length(r),ncol = 1) # unit vector

  a = drop(t(e)%*% invS %*% e)
  b = drop(t(e)%*% invS %*% r)
  c = drop(t(r)%*% invS %*% r)
  d = drop(a*c - b^2)
  
  
  yy= seq(from = min_r, to = max_r,length.out = N+1)
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





EfficientFrontier_constr = function(r,S,full=FALSE,plot=FALSE, N=100, no_short_sales, max_r=NA, min_r =NA){
  # INPUT
  # r: expected returns of the assets
  # S: covariance matrix of the asset returns
  # full: computes only upper section of frontier if FALSE
  # plot: if TRUE the frontier is plotted
  # N: how many expected returns to take into consideration to plot frontier

  library(quadprog)
  
  if (is.na(max_r)){
    max_r = max(r)
  }
  if (is.na(min_r)){
    min_r = min(r)
  }
  
  # number of assets
  n= length(r)

  D = 2*S       #times 2 because there is a 1/2 in the implicit formulation
  d = rep(0,n)  # zeros
  
  
  # Constraint on returns
  A = t(r)
  
  b = 0.0 # expected return that will be iterated to compute frontier
  
  # Constraint on weights: sum(w)=1
  A = rbind(A,rep(1,n))
  b= rbind(b,1)
  
  
  # Constraint on short selling ony on given assets
  for (i in no_short_sales){
    A_i = rep(0,n)
    A_i[i]=1
    A = rbind(A,A_i)
    b = rbind(b,0)
  }

  # yy = set of returns
  # xx = set of corresponding volatility
  yy= seq(from = min(r), to = max_r,length.out = N+1)
  yy[1] = yy[1]+ 1e-5
  yy[N+1] = yy[N+1] -1e-5
  
  b[1]= yy[1]
  # print(yy)
  
  xx = rep(0,length(yy))
  for (i in 1:(N+1)) {
    # print(i)
    # print(yy[i])
    b[1] = yy[i]
    # print(b)
    sol = solve.QP(Dmat = D, dvec = (d), Amat = t(A), bvec = t(b), meq = 2)
    xx[i] = sqrt(sol$value)
  }

  if (!full){
    min_sigma = min(xx)
    # print(min_sigma)
    idx = which(yy>=yy[which(xx==min_sigma)])
    # print(idx)
    xx= xx[idx]
    yy= yy[idx]
  }

  # if(plot){
  #   min_y = min(c(yy,r))
  #   max_y = max(c(yy,r))
  #   plot(xx,yy,type ='l', ylim = c(min_y,max_y))
  #   points(sqrt(diag(S)),r,col='blue',pch='+')
  #   #legend(c("Efficient Frontier","Single Assets"))
  # }
  res = list(sigma = xx, expected_return=yy)
  return(res)
}
