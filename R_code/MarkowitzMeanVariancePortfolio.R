###############################################################
############## MARKOVITZ MEAN VARIANCE PORTFOLIO ##############
###############################################################

<<<<<<< HEAD
<<<<<<< HEAD
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
  
  yy= seq(from = -0.1, to = 0.2,length.out = N+1)
  xx = sqrt( (a*yy^2 - 2*b*yy + c)/d)
  
  if (!full){
    min_sigma = min(xx)
    idx = which(yy>=yy[which(xx==min_sigma)])
    # print(idx)
    xx= xx[idx]
    yy= yy[idx]
=======
=======
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
# Simple wrapper for frontier
EfficientFrontier=function(r,S, no_short_sales=NA, N=100,full=FALSE,plot=FALSE, max_r=NA){
  
  if (is.na(no_short_sales[1])){
    return(EfficientFrontier_unconstr(r=r,S=S,full=full,plot=plot, N=N))
<<<<<<< HEAD
>>>>>>> Added constr and unconstr optimal allocation, added plots.
=======
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
  }
  else{
    return(EfficientFrontier_constr(r=r,S=S,full=full,plot=plot, N=N, no_short_sales = no_short_sales))
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



<<<<<<< HEAD
<<<<<<< HEAD

OptimalAllocation= function(r,S, expected_return = NA, sd = NA){
  
  if(is.na(expected_return) & is.na(sd)){
=======
OptimalAllocation_unconstr= function(r,S, target_return = NA, sd = NA){
  
  if(is.na(target_return) & is.na(sd)){
>>>>>>> Added constr and unconstr optimal allocation, added plots.
=======
OptimalAllocation_unconstr= function(r,S, target_return = NA, sd = NA){
  
  if(is.na(target_return) & is.na(sd)){
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
    stop("Please specify expected return or standard deviation.")
  }
  
  invS = solve(S)
<<<<<<< HEAD
<<<<<<< HEAD
  e = matrix(rep(1,length(r)),nrow = length(r),ncol = 1) # unit vector 
=======
  e = matrix(rep(1,length(r)),nrow = length(r),ncol = 1) # unit vector
>>>>>>> Added constr and unconstr optimal allocation, added plots.
=======
  e = matrix(rep(1,length(r)),nrow = length(r),ncol = 1) # unit vector
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
  
  a = drop(t(e)%*% invS %*% e)
  b = drop(t(e)%*% invS %*% r)
  c = drop(t(r)%*% invS %*% r)
  d = drop(a*c - b^2)
  
  
<<<<<<< HEAD
<<<<<<< HEAD
  if (!is.na(expected_return)){
    mu = drop((a*expected_return - b)/d)
    lambda = drop( (d - a*b*expected_return + b^2)/(a*d))
=======
  if (!is.na(target_return)){
    mu = drop((a*target_return - b)/d)
    lambda = drop( (d - a*b*target_return + b^2)/(a*d))
>>>>>>> Added constr and unconstr optimal allocation, added plots.
=======
  if (!is.na(target_return)){
    mu = drop((a*target_return - b)/d)
    lambda = drop( (d - a*b*target_return + b^2)/(a*d))
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
    
    w_opt = invS %*% (e*lambda + mu*r)
    
    sigma_opt = sqrt(t(w_opt)%*%S%*%w_opt)
<<<<<<< HEAD
<<<<<<< HEAD
    sigma = sqrt( (a*expected_return^2 - 2*b*expected_return + c)/d)
=======
    sigma = sqrt( (a*target_return^2 - 2*b*target_return + c)/d)
>>>>>>> Added constr and unconstr optimal allocation, added plots.
=======
    sigma = sqrt( (a*target_return^2 - 2*b*target_return + c)/d)
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
    
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
<<<<<<< HEAD
<<<<<<< HEAD
   
=======
    
>>>>>>> Added constr and unconstr optimal allocation, added plots.
=======
    
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
    
    # test on the results:
    if (abs(sd - sigma_opt) > 1e-6){
      stop("Computation of minimum standard deviation yields a different result.")
    }
    print(paste("Maximum expected return for given standard deviation is",expected_return_opt))
  }
  
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
  return(w_opt)
}



OptimalAllocation_constr= function(r,S, target_return = NA, no_short_sales){
  
  if(is.na(target_return)){
    stop("Please specify expected return.")
  }
  if(target_return>max(r)){
    warning("Cannot produce such a high return without short-selling.")
    target_return=max(r)
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
<<<<<<< HEAD
>>>>>>> Added constr and unconstr optimal allocation, added plots.
=======
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
  return(w_opt)
}



<<<<<<< HEAD
<<<<<<< HEAD
EfficientFrontier_constr = function(r,S,full=TRUE,plot=TRUE, N=100){
  # INPUT 
=======
=======
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349


EfficientFrontier_unconstr = function(r,S,full=FALSE,plot=FALSE, N=100, max_r=NA){
  # INPUT
<<<<<<< HEAD
>>>>>>> Added constr and unconstr optimal allocation, added plots.
=======
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
  # r: expected returns of the assets
  # S: covariance matrix of the asset returns
  # full: computes only upper section of frontier if FALSE
  # plot: if TRUE the frontier is plotted
  # N: how many expected returns to take into consideration to plot frontier
<<<<<<< HEAD
<<<<<<< HEAD
  
  library(quadprog)
  
  
  
  D = 2*S
  d = matrix(rep(0,length(r)),ncol = 1)
  A = rbind(t(r),rep(1,length(r)),
            diag(length(r)))
  
  
  sol = solve.QP(Dmat = D, dvec = (d), Amat = t(A), bvec = t(b), meq = 2)
  
  yy= seq(from = min(r), to = max(r),length.out = N+1)
  b= c(yy[1], 1, rep(0,length(r)))
=======
=======
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
  if (is.na(max_r)){
    max_r = max(r)
  }

  invS = solve(S)
  e = matrix(rep(1,length(r)),nrow = length(r),ncol = 1) # unit vector

  a = drop(t(e)%*% invS %*% e)
  b = drop(t(e)%*% invS %*% r)
  c = drop(t(r)%*% invS %*% r)
  d = drop(a*c - b^2)
  
  
  yy= seq(from = -0.1, to = max_r,length.out = N+1)
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





EfficientFrontier_constr = function(r,S,full=FALSE,plot=FALSE, N=100, no_short_sales, max_r=NA){
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


  # yy= set of returns
  # xx =set of corresponding volatility
  yy= seq(from = min(r), to = max_r,length.out = N+1)
  yy[1] = yy[1]+ 1e-5
  
  b[1]= yy[1]

  
<<<<<<< HEAD
>>>>>>> Added constr and unconstr optimal allocation, added plots.
  xx = rep(0,length(yy))
  for (i in 1:(N+1)) {
    b = c(yy[i], 1, rep(0,length(r)))
    sol = solve.QP(Dmat = D, dvec = (d), Amat = t(A), bvec = t(b), meq = 2)
    xx[i] = sqrt(sol$value)
  }
<<<<<<< HEAD
  
  # if (!full){
  #   min_sigma = min(xx)
  #   idx = which(yy>=yy[which(xx==min_sigma)])
  #   # print(idx)
  #   xx= xx[idx]
  #   yy= yy[idx]
  # }
  
=======
=======
  xx = rep(0,length(yy))
  for (i in 1:(N+1)) {
    #print(yy[i])
    b[1] = yy[i]
    sol = solve.QP(Dmat = D, dvec = (d), Amat = t(A), bvec = t(b), meq = 2)
    xx[i] = sqrt(sol$value)
  }
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349

  if (!full){
    min_sigma = min(xx)
    # print(min_sigma)
    idx = which(yy>=yy[which(xx==min_sigma)])
    # print(idx)
    xx= xx[idx]
    yy= yy[idx]
  }

<<<<<<< HEAD
>>>>>>> Added constr and unconstr optimal allocation, added plots.
=======
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
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
<<<<<<< HEAD

=======
>>>>>>> 9387df364c7056a942a308a836cc64b4c5b77349
