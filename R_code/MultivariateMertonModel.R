############################################################
################ Multivariate Merton Model #################
############################################################
# This library computes the transition density of a Multivariate 
# Merton model in different cases and implements the loglikelihood function. 
# Explicit implementations are needed to make the code run much faster,
# that's why there are explicit functions for each number of asset up to 4


library(mvtnorm)





MultivariateMertonPdf = function(x, dt, mu, S, theta, delta, lambda, theta_z, delta_z, lambda_z, alpha){
  # Computes the density of a multivariate merton model returns with idiosyncratic and common jumps
  # NOTE: all vectors should be vertical [n*1]
  # ASSUMPTION: in dt time we can only have 0 or 1 jumps in each jump process, so lambda*dt<=1
  #
  # INPUT
  # x:      vector representing at which point to compute the density [vector of n]
  # mu:     drift of the continuos part [vector of n]
  # S:      covariance of the continuous part [matrix n*n]
  # theta:  means of the idiosyncratic jump intensity [vector of n]
  # delta:  variances of the idiosyncratic jump intensity [vector of n]
  # lambda:  poisson parameters of the idiosyncratic jump part [vector of n]
  # theta_z: mean of common jump intensity 
  # delta_z: variance of common jump intensity
  # lambda_z: poisson parameter of common jump part
  # alpha:  vector of coefficient that multiply the common jump effect for each component
  
  # check on lambdas: 
  require(binaryLogic)
  
  if(sum(lambda*dt>=1) || lambda_z*dt >=1){
    stop("Error: lambda*dt should be lower than 1 (ideally closer to 0).")
  }
  
  n = length(mu)
  
  cov_z = alpha%*%t(alpha)*delta_z
  mu_z = theta_z*alpha
  
  pdf = 0
  
  for(k in 0:(2^(n+1)-1)){
    mean_x = mu*dt
    cov_x = S*sqrt(dt)
    prob = 1
    
    k_bin = as.binary(k,n = n+1,littleEndian=T)
    
    # building probability and conditional density
    for (i in 1:(n+1)){
      
      if ( i<n+1){
        if ( k_bin[i]){
          prob = prob*lambda[i]*dt
          mean_x[i] = mean_x[i] + theta[i]
          cov_x[i,i] = cov_x[i,i] + delta[i]^2
        }
        else{
          prob = prob*(1-lambda[i]*dt)
        }
      }
      
      else{
        if ( k_bin[i]){
          prob = prob*lambda_z*dt
          mean_x = mean_x + mu_z
          cov_x = cov_x + cov_z
        }
        else
          prob = prob*(1-lambda_z*dt)
      }
    }
    
    # adding each term
    # print(length(x))
    partial_pdf = dmvnorm(x,mean = mean_x, sigma = cov_x)
    pdf = pdf + prob*partial_pdf
  }
  
  return(pdf)
}



negloglik = function(params, x, dt, n) {
  # 
  # x is a matrix [Npoints * n] of all the points for which we compute the likelihood
  # 
  
  
  ## add check on inputs
  
  # reconstruction of parameters:
  idx =1
  mu = params[idx:(idx+n-1)]
  idx = idx+n 
  
  
  S = matrix(rep(0,n*n), ncol = n)
  i=1
  j=1
  for(k in 1:(n*(n+1)/2)){
    S[i,j] = params[idx+k-1]
    S[j,i] =  S[i,j]
    j=j+1
    if(j == n+1){
      i=i+1
      j=i
    }
  }
  idx = idx + n*(n+1)/2
  
  theta = params[idx:(idx+n-1)]
  idx = idx+n
  
  delta = params[idx:(idx+n-1)]
  idx = idx+n
  
  lambda = params[idx:(idx+n-1)]
  idx = idx+n
  
  theta_z = params[idx]
  idx = idx+1
  
  delta_z = params[idx]
  idx = idx+1
  
  lambda_z = params[idx]
  idx = idx+1
  
  alpha = params[idx:(idx+n-1)]
  idx = idx+n
  
  # print(mu)
  # print(S)
  # print(theta)
  # print(delta)
  # print(lambda)
  # print(alpha)
  if( (idx-1)!=length(params))
    stop("Error in parameter reconstruction: number of parameters is wrong.")
  
  
  # computing pdf on each point and adding
  partial = 0
  for(i in 1:dim(x)[1]){
    pdf = MultivariateMertonPdf(x[i,], dt, mu, S, theta, delta, lambda, theta_z, delta_z, lambda_z, alpha)
    # cat("\npdf:")
    # print(pdf)
    partial = partial + log(pdf)
  }
  
  # last check on results
  nll = -(partial)
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }
  return(nll)
}








#############################################################################
################################### no common ###############################
#############################################################################

MultivariateMertonPdf_nocommon = function(x, dt, mu, S, theta, delta, lambda){
  # Computes the density of a multivariate merton model returns with idiosyncratic and common jumps
  # NOTE: all vectors should be vertical [n*1]
  # ASSUMPTION: in dt time we can only have 0 or 1 jumps in each jump process, so lambda*dt<=1
  #
  # INPUT
  # x:      vector representing at which point to compute the density [vector of n]
  # mu:     drift of the continuos part [vector of n]
  # S:      covariance of the continuous part [matrix n*n]
  # theta:  means of the idiosyncratic jump intensity [vector of n]
  # delta:  variances of the idiosyncratic jump intensity [vector of n]
  # lambda:  poisson parameters of the idiosyncratic jump part [vector of n]
  # theta_z: mean of common jump intensity 
  # delta_z: variance of common jump intensity
  # lambda_z: poisson parameter of common jump part
  # alpha:  vector of coefficient that multiply the common jump effect for each component
  
  # check on lambdas: 
  require(binaryLogic)
  
  if(sum(lambda*dt>=1)){
    stop("Error: lambda*dt should be lower than 1 (ideally closer to 0).")
  }
  
  n = length(mu)
  
  pdf = 0
  
  for(k in 0:(2^(n)-1)){
    mean_x = mu*dt
    cov_x = S*sqrt(dt)
    prob = 1
    
    k_bin = as.binary(k,n = n+1,littleEndian=T)
    
    # building probability and conditional density
    for (i in 1:(n)){
        if ( k_bin[i]){
          prob = prob*lambda[i]*dt
          mean_x[i] = mean_x[i] + theta[i]
          cov_x[i,i] = cov_x[i,i] + delta[i]^2
        }
        else{
          prob = prob*(1-lambda[i]*dt)
        }
    }
    
    # adding each term
    # print(length(x))
    
    partial_pdf = dmvnorm(x,mean = mean_x, sigma = cov_x)
    pdf = pdf + prob*partial_pdf
    # print(prob)
    # print(partial_pdf)
    # print(pdf)
  }

  return(pdf)
}

negloglik_nocommon = function(params, x, dt, n) {
  # 
  # x is a matrix [Npoints * n] of all the points for which we compute the likelihood
  # 

  ## add check on inputs
  
  # reconstruction of parameters:
  idx =1
  mu = params[idx:(idx+n-1)]
  idx = idx+n 


  sigma=diag(params[idx:(idx+n-1)])
  idx = idx + n
  
  corr = matrix(rep(0,n*n), ncol = n)
  k=1
  for(i in 1:n)
    for(j in i:n){
      if (j!=i){
        corr[i,j]=params[idx+k-1]
        corr[j,i]=params[idx+k-1]
        k=k+1
      }
      else
        corr[i,j]=1
    }
  idx = idx + n*(n-1)/2
  
  S = sigma %*% corr %*% sigma
  
  # print(S)
  
  theta = params[idx:(idx+n-1)]
  idx = idx+n
  
  delta = params[idx:(idx+n-1)]
  idx = idx+n
  
  lambda = params[idx:(idx+n-1)]
  idx = idx+n
  
  # print(mu)
  # print(S)
  # print(theta)
  # print(delta)
  # print(lambda)
  # print(alpha)
  if( (idx-1)!=length(params))
    stop("Error in parameter reconstruction: number of parameters is wrong.")
  
  
  # computing pdf on each point and adding
  partial = 0
  for(i in 1:dim(x)[1]){
    pdf = MultivariateMertonPdf_nocommon(x[i,], dt, mu, S, theta, delta, lambda)
    # cat("\npdf:")
    # print(pdf)
    partial = partial + log(pdf)
  }
  
  # last check on results
  nll = -(partial)
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }
  return(nll)
}


#############################################################################
############################### one asset  ##################################
#############################################################################


MultivariateMertonPdf_1asset_nocommon = function(x, dt, mu, S, theta, delta, lambda, a=T, b=T){
  # Computes the density of a multivariate merton model returns with idiosyncratic and common jumps
  # NOTE: all vectors should be vertical [n*1]
  # ASSUMPTION: in dt time we can only have 0 or 1 jumps in each jump process, so lambda*dt<=1
  #
  # INPUT
  # x:      vector representing at which point to compute the density [vector of n]
  # mu:     drift of the continuos part [vector of n]
  # S:      covariance of the continuous part [matrix n*n]
  # theta:  means of the idiosyncratic jump intensity [vector of n]
  # delta:  variances of the idiosyncratic jump intensity [vector of n]
  # lambda:  poisson parameters of the idiosyncratic jump part [vector of n]

  # check on lambdas: 
  ldt =lambda*dt

  if(ldt>=1){
    stop("Error: lambda*dt should be lower than 1 (ideally close to 0).")
  }
  mu_adj= mu -lambda*theta - S/2#
  pdf=a*(1-ldt)*dnorm(x, mean = mu_adj*dt, sd = sqrt(S*dt))+
      b*ldt *dnorm(x, mean = mu_adj*dt + log(1+theta)-S/2, sd = sqrt(S*dt + delta^2))
  
  #print(ldt)
  return(pdf)
}


negloglik_1asset_nocommon= function(params, x, dt, n) {
  # 
  # x is a matrix [Npoints * n] of all the points for which we compute the likelihood
  # 
  ## add check on inputs
  
  # reconstruction of parameters:
  mu=params[1]
  S = params[2]^2
  theta = params[3]
  delta = params[4]^2
  lambda = params[5]
  
  
  # print(mu)
  # print(S)
  # print(theta)
  # print(delta)
  # print(lambda)
  # print(alpha)

  
  
  # computing pdf on each point and adding
  partial = MultivariateMertonPdf_1asset_nocommon(x, dt, mu, S, theta, delta, lambda)
  nll = -sum(log(partial))
  
  # last check on result
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }
  
  #print("ONEASSETDENSITY")
  return(nll)
}


#############################################################################
####################### two assets + 1 common jump ##########################
#############################################################################


MultivariateMertonPdf_2assets = function(x, dt, mu, S, theta, delta, lambda, theta_z, delta_z, lambda_z, alpha){
  # Computes the density of a multivariate merton model returns with idiosyncratic and common jumps
  # NOTE: all vectors should be vertical [n*1]
  # ASSUMPTION: in dt time we can only have 0 or 1 jumps in each jump process, so lambda*dt<=1
  #
  # INPUT
  # x:      vector representing at which point to compute the density [vector of n]
  # mu:     drift of the continuos part [vector of n]
  # S:      covariance of the continuous part [matrix n*n]
  # theta:  means of the idiosyncratic jump intensity [vector of n]
  # delta:  variances of the idiosyncratic jump intensity [vector of n]
  # lambda:  poisson parameters of the idiosyncratic jump part [vector of n]
  # theta_z: mean of common jump intensity 
  # delta_z: variance of common jump intensity
  # lambda_z: poisson parameter of common jump part
  # alpha:  vector of coefficient that multiply the common jump effect for each component
  
  # check on lambdas: 
  ldt =lambda*dt
  ldt_z = lambda_z*dt
  
  if(sum(ldt>=1) || ldt_z >=1){
    stop("Error: lambda*dt should be lower than 1 (ideally closer to 0).")
  }
  
  n = length(mu)
  
  cov_z = alpha%*%t(alpha)*delta_z
  mean_z = theta_z*alpha
  
  pdf = 0
  
  mean_y1 = c(theta[1],0)
  cov_y1 = matrix(c(delta[1]^2,0,0,0), 2,2)
  
  mean_y2 = c(0,theta[2])
  cov_y2 = matrix(c(0,0,0,delta[2]^2),2,2)
  
  mean_x = mu*dt
  cov_x = S*(dt)
  

  #000
  pdf= pdf + (1-ldt[1])*(1-ldt[2])*(1-ldt_z)*dmvnorm(x, mean = mean_x, sigma = cov_x)
  #001
  pdf= pdf + (ldt[1])*(1-ldt[2])*(1-ldt_z)*dmvnorm(x, mean = mean_x+mean_y1, sigma = cov_x+cov_y1)
  #010
  pdf= pdf + (1-ldt[1])*(ldt[2])*(1-ldt_z)*dmvnorm(x, mean = mean_x+mean_y2, sigma = cov_x+cov_y2)
  #011
  pdf= pdf + (ldt[1])*(ldt[2])*(1-ldt_z)*dmvnorm(x, mean = mean_x+mean_y1+mean_y2, sigma = cov_x+cov_y1+cov_y2)
  #100
  pdf= pdf + (1-ldt[1])*(1-ldt[2])*(ldt_z)*dmvnorm(x, mean = mean_x+mean_z, sigma = cov_x+cov_z)
  #101
  pdf= pdf + (ldt[1])*(1-ldt[2])*(ldt_z)*dmvnorm(x, mean = mean_x+mean_y1+mean_z, sigma = cov_x+cov_y1+cov_z)
  #110
  pdf= pdf + (1-ldt[1])*(ldt[2])*(ldt_z)*dmvnorm(x, mean = mean_x+mean_y2+mean_z, sigma = cov_x+cov_y2+cov_z)
  #111
  pdf= pdf + (ldt[1])*(ldt[2])*(ldt_z)*dmvnorm(x, mean = mean_x+mean_y1+mean_y2+mean_z, sigma = cov_x+cov_y1+cov_y2+cov_z)
  
  
  return(pdf)
}


negloglik_2assets= function(params, x, dt, n) {
  # 
  # x is a matrix [Npoints * n] of all the points for which we compute the likelihood
  # 
  ## add check on inputs
  
  # reconstruction of parameters:
  idx =1
  mu = params[idx:(idx+n-1)]
  idx = idx+n 
  
  
  sigma=diag(params[idx:(idx+n-1)])
  idx = idx + n
  
  corr = matrix(rep(0,n*n), ncol = n)
  k=1
  for(i in 1:n)
    for(j in i:n){
      if (j!=i){
        corr[i,j]=params[idx+k-1]
        corr[j,i]=params[idx+k-1]
        k=k+1
      }
      else
        corr[i,j]=1
    }
  idx = idx + n*(n-1)/2
  
  S = sigma %*% corr %*% sigma
  
  theta = params[idx:(idx+n-1)]
  idx = idx+n
  
  delta = params[idx:(idx+n-1)]
  idx = idx+n
  
  lambda = params[idx:(idx+n-1)]
  idx = idx+n
  
  theta_z = params[idx]
  idx = idx+1
  
  delta_z = params[idx]
  idx = idx+1
  
  lambda_z = params[idx]
  idx = idx+1
  
  alpha = params[idx:(idx+n-1)]
  idx = idx+n
  
  # print(mu)
  # print(S)
  # print(theta)
  # print(delta)
  # print(lambda)
  # print(alpha)
  if( (idx-1)!=length(params))
    stop("Error in parameter reconstruction: number of parameters is wrong.")
  
  
  # computing pdf on each point and adding
  partial = MultivariateMertonPdf_2assets(x, dt, mu, S, theta, delta, lambda, theta_z, delta_z, lambda_z, alpha)
  nll = -sum(log(partial))
  
  # last check on result
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }
  return(nll)
}


#############################################################################
####################### two assets no common jump ###########################
#############################################################################


MultivariateMertonPdf_2assets_nocommon = function(x, dt, mu, S, theta, delta, lambda){
  # Computes the density of a multivariate merton model returns with idiosyncratic and common jumps
  # NOTE: all vectors should be vertical [n*1]
  # ASSUMPTION: in dt time we can only have 0 or 1 jumps in each jump process, so lambda*dt<=1
  #
  # INPUT
  # x:      vector representing at which point to compute the density [vector of n]
  # mu:     drift of the continuos part [vector of n]
  # S:      covariance of the continuous part [matrix n*n]
  # theta:  means of the idiosyncratic jump intensity [vector of n]
  # delta:  variances of the idiosyncratic jump intensity [vector of n]
  # lambda:  poisson parameters of the idiosyncratic jump part [vector of n]
  # theta_z: mean of common jump intensity 
  # delta_z: variance of common jump intensity
  # lambda_z: poisson parameter of common jump part
  # alpha:  vector of coefficient that multiply the common jump effect for each component
  
  # check on lambdas: 
  ldt =lambda*dt
  
  if(sum(ldt>=1)){
    stop("Error: lambda*dt should be lower than 1 (ideally closer to 0).")
  }
  
  n = length(mu)
  
  pdf = 0
  
  mean_y1 = c(log(1+theta[1])-0.5 * delta[1],0)
  cov_y1 = matrix(c(delta[1]^2,0,0,0), 2,2)
  
  mean_y2 = c(0,log(1+theta[2])-0.5 * delta[2])
  cov_y2 = matrix(c(0,0,0,delta[2]^2),2,2)
  
  mean_x = mu*dt - ldt * theta - diag(S)/2
  cov_x = S*(dt)
  
  
  #000
  pdf= pdf + (1-ldt[1])*(1-ldt[2])*dmvnorm(x, mean = mean_x, sigma = cov_x)
  #001
  pdf= pdf + (ldt[1])*(1-ldt[2])*dmvnorm(x, mean = mean_x+mean_y1, sigma = cov_x+cov_y1)
  #010
  pdf= pdf + (1-ldt[1])*(ldt[2])*dmvnorm(x, mean = mean_x+mean_y2, sigma = cov_x+cov_y2)
  #011
  pdf= pdf + (ldt[1])*(ldt[2])*dmvnorm(x, mean = mean_x+mean_y1+mean_y2, sigma = cov_x+cov_y1+cov_y2)
  
  return(pdf)
}


negloglik_2assets_nocommon= function(params, x, dt, n) {
  # 
  # x is a matrix [Npoints * n] of all the points for which we compute the likelihood
  # 
  ## add check on inputs
  # reconstruction of parameters:
  idx =1
  mu = params[idx:(idx+n-1)]
  idx = idx+n 
  
  
  sigma=diag(params[idx:(idx+n-1)])
  idx = idx + n
  
  corr = matrix(rep(0,n*n), ncol = n)
  k=1
  for(i in 1:n)
    for(j in i:n){
      if (j!=i){
        corr[i,j]=params[idx+k-1]
        corr[j,i]=params[idx+k-1]
        k=k+1
      }
      else
        corr[i,j]=1
    }
  idx = idx + n*(n-1)/2
  
  S = sigma %*% corr %*% sigma
  
  
  theta = params[idx:(idx+n-1)]
  idx = idx+n
  
  delta = params[idx:(idx+n-1)]
  idx = idx+n
  
  lambda = params[idx:(idx+n-1)]
  idx = idx+n
  
  # print(mu)
  # print(S)
  # print(theta)
  # print(delta)
  # print(lambda)

  if( (idx-1)!=length(params))
    stop("Error in parameter reconstruction: number of parameters is wrong.")
  
  
  # computing pdf on each point and adding
  partial = MultivariateMertonPdf_2assets_nocommon(x, dt, mu, S, theta, delta, lambda)
  nll = -sum(log(partial))
  
  # last check on result
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }
  return(nll)
}


#############################################################################
####################### three assets no common jump #########################
#############################################################################


MultivariateMertonPdf_3assets_nocommon = function(x, dt, mu, S, theta, delta, lambda){
  # Computes the density of a multivariate merton model returns with idiosyncratic and common jumps
  # NOTE: all vectors should be vertical [n*1]
  # ASSUMPTION: in dt time we can only have 0 or 1 jumps in each jump process, so lambda*dt<=1
  #
  # INPUT
  # x:      vector representing at which point to compute the density [vector of n]
  # mu:     drift of the continuos part [vector of n]
  # S:      covariance of the continuous part [matrix n*n]
  # theta:  means of the idiosyncratic jump intensity [vector of n]
  # delta:  variances of the idiosyncratic jump intensity [vector of n]
  # lambda:  poisson parameters of the idiosyncratic jump part [vector of n]
  # theta_z: mean of common jump intensity 
  # delta_z: variance of common jump intensity
  # lambda_z: poisson parameter of common jump part
  # alpha:  vector of coefficient that multiply the common jump effect for each component
  
  # check on lambdas: 
  ldt =lambda*dt
  
  if(sum(ldt>=1)){
    stop("Error: lambda*dt should be lower than 1 (ideally closer to 0).")
  }
  
  n = length(mu)
  
  pdf = 0
  
  mean_y1 = c(theta[1],0,0)
  cov_y1 = matrix(c(delta[1]^2,0,0,0,0,0,0,0,0), 3,3)
  
  mean_y2 = c(0,theta[2],0)
  cov_y2 = matrix(rep(0,3*3),3,3)
  cov_y2[2,2] = delta[2]^2
  
  mean_y3 = c(0,0,theta[3])
  cov_y3 = matrix(rep(0,3*3),3,3)
  cov_y3[3,3] = delta[3]^2
  
  mean_x = mu*dt
  cov_x = S*(dt)
  
  
  #000
  pdf= pdf + (1-ldt[1])*(1-ldt[2])*(1-ldt[3])*dmvnorm(x, mean = mean_x, sigma = cov_x)
  #001
  pdf= pdf + (ldt[1])*(1-ldt[2])*(1-ldt[3])*dmvnorm(x, mean = mean_x+mean_y1, sigma = cov_x+cov_y1)
  #010
  pdf= pdf + (1-ldt[1])*(ldt[2])*(1-ldt[3])*dmvnorm(x, mean = mean_x+mean_y2, sigma = cov_x+cov_y2)
  #011
  pdf= pdf + (ldt[1])*(ldt[2])*(1-ldt[3])*dmvnorm(x, mean = mean_x+mean_y1+mean_y2, sigma = cov_x+cov_y1+cov_y2)
  #100
  pdf= pdf + (1-ldt[1])*(1-ldt[2])*(ldt[3])*dmvnorm(x, mean = mean_x+mean_y3, sigma = cov_x+cov_y3)
  #101
  pdf= pdf + (ldt[1])*(1-ldt[2])*(ldt[3])*dmvnorm(x, mean = mean_x+mean_y1+mean_y3, sigma = cov_x+cov_y1+cov_y3)
  #110
  pdf= pdf + (1-ldt[1])*(ldt[2])*(ldt[3])*dmvnorm(x, mean = mean_x+mean_y2+mean_y3, sigma = cov_x+cov_y2+cov_y3)
  #111
  pdf= pdf + (ldt[1])*(ldt[2])*(ldt[3])*dmvnorm(x, mean = mean_x+mean_y1+mean_y2+mean_y3, sigma = cov_x+cov_y1+cov_y2+cov_y3)
  
  return(pdf)
}


negloglik_3assets_nocommon= function(params, x, dt, n) {
  # 
  # x is a matrix [Npoints * n] of all the points for which we compute the likelihood
  # 
  ## add check on inputs
  
  # reconstruction of parameters:
  idx =1
  mu = params[idx:(idx+n-1)]
  idx = idx+n 
  
  
  sigma=diag(params[idx:(idx+n-1)])
  idx = idx + n
  
  corr = matrix(rep(0,n*n), ncol = n)
  k=1
  for(i in 1:n)
    for(j in i:n){
      if (j!=i){
        corr[i,j]=params[idx+k-1]
        corr[j,i]=params[idx+k-1]
        k=k+1
      }
      else
        corr[i,j]=1
    }
  idx = idx + n*(n-1)/2
  
  S = sigma %*% corr %*% sigma
  
  theta = params[idx:(idx+n-1)]
  idx = idx+n
  
  delta = params[idx:(idx+n-1)]
  idx = idx+n
  
  lambda = params[idx:(idx+n-1)]
  idx = idx+n
  
  # print(mu)
  # print(S)
  # print(theta)
  # print(delta)
  # print(lambda)
  # print(alpha)
  if( (idx-1)!=length(params))
    stop("Error in parameter reconstruction: number of parameters is wrong.")
  
  
  # computing pdf on each point and adding
  partial = MultivariateMertonPdf_3assets_nocommon(x, dt, mu, S, theta, delta, lambda)
  nll = -sum(log(partial))
  
  # last check on result
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }
  return(nll)
}

#############################################################################
####################### four assets no common jump ##########################
#############################################################################


MultivariateMertonPdf_4assets_nocommon = function(x, dt, mu, S, theta, delta, lambda){
  # Computes the density of a multivariate merton model returns with idiosyncratic and common jumps
  # NOTE: all vectors should be vertical [n*1]
  # ASSUMPTION: in dt time we can only have 0 or 1 jumps in each jump process, so lambda*dt<=1
  #
  # INPUT
  # x:      vector representing at which point to compute the density [vector of n]
  # mu:     drift of the continuos part [vector of n]
  # S:      covariance of the continuous part [matrix n*n]
  # theta:  means of the idiosyncratic jump intensity [vector of n]
  # delta:  variances of the idiosyncratic jump intensity [vector of n]
  # lambda:  poisson parameters of the idiosyncratic jump part [vector of n]
  # theta_z: mean of common jump intensity 
  # delta_z: variance of common jump intensity
  # lambda_z: poisson parameter of common jump part
  # alpha:  vector of coefficient that multiply the common jump effect for each component
  
  # check on lambdas: 
  ldt =lambda*dt
  
  if(sum(ldt>=1)){
    stop("Error: lambda*dt should be lower than 1 (ideally closer to 0).")
  }
  
  n = length(mu)
  
  pdf = 0
  
  mean_y1 = rep(0,n)
  mean_y1[1] = theta[1]
  cov_y1 = matrix(rep(0,n*n),n,n)
  cov_y1[1,1] = delta[1]^2

  mean_y2 = rep(0,n)
  mean_y2[2] = theta[2]
  cov_y2 = matrix(rep(0,n*n),n,n)
  cov_y2[2,2] = delta[2]^2
  
  mean_y3 = rep(0,n)
  mean_y3[3] = theta[3]
  cov_y3 = matrix(rep(0,n*n),n,n)
  cov_y3[3,3] = delta[3]^2
  
  mean_y4 = rep(0,n)
  mean_y4[4] = theta[4]
  cov_y4 = matrix(rep(0,n*n),n,n)
  cov_y4[4,4] = delta[4]^2
  
  mean_x = (mu - S[c(1,6,11,16)]*0.5 ) *dt #- lambda*theta*dt #FIXME ldt o non ldt?
  cov_x = S*(dt) # FIXED: check sqrt(dt) or dt -->FIXED dt because it is a covariance not a volatility
  
  #0000
  pdf= pdf + (1-ldt[1])*(1-ldt[2])*(1-ldt[3])*(1-ldt[4])*dmvnorm(x, mean = mean_x, sigma = cov_x)
  #0001
  pdf= pdf + (ldt[1])*(1-ldt[2])*(1-ldt[3])*(1-ldt[4])*dmvnorm(x, mean = mean_x+mean_y1, sigma = cov_x+cov_y1)
  #0010
  pdf= pdf + (1-ldt[1])*(ldt[2])*(1-ldt[3])*(1-ldt[4])*dmvnorm(x, mean = mean_x+mean_y2, sigma = cov_x+cov_y2)
  #0011
  pdf= pdf + (ldt[1])*(ldt[2])*(1-ldt[3])*(1-ldt[4])*dmvnorm(x, mean = mean_x+mean_y1+mean_y2, sigma = cov_x+cov_y1+cov_y2)
  #0100
  pdf= pdf + (1-ldt[1])*(1-ldt[2])*(ldt[3])*(1-ldt[4])*dmvnorm(x, mean = mean_x+mean_y3, sigma = cov_x+cov_y3)
  #0101
  pdf= pdf + (ldt[1])*(1-ldt[2])*(ldt[3])*(1-ldt[4])*dmvnorm(x, mean = mean_x+mean_y1+mean_y3, sigma = cov_x+cov_y1+cov_y3)
  #0110
  pdf= pdf + (1-ldt[1])*(ldt[2])*(ldt[3])*(1-ldt[4])*dmvnorm(x, mean = mean_x+mean_y2+mean_y3, sigma = cov_x+cov_y2+cov_y3)
  #0111
  pdf= pdf + (ldt[1])*(ldt[2])*(ldt[3])*(1-ldt[4])*dmvnorm(x, mean = mean_x+mean_y1+mean_y2+mean_y3, sigma = cov_x+cov_y1+cov_y2+cov_y3)
  #1000
  pdf= pdf + (1-ldt[1])*(1-ldt[2])*(1-ldt[3])*(ldt[4])*dmvnorm(x, mean = mean_x+mean_y4, sigma = cov_x+cov_y4)
  #1001
  pdf= pdf + (ldt[1])*(1-ldt[2])*(1-ldt[3])*(ldt[4])*dmvnorm(x, mean = mean_x+mean_y1+mean_y4, sigma = cov_x+cov_y1+cov_y4)
  #1010
  pdf= pdf + (1-ldt[1])*(ldt[2])*(1-ldt[3])*(ldt[4])*dmvnorm(x, mean = mean_x+mean_y2+mean_y4, sigma = cov_x+cov_y2+cov_y4)
  #1011
  pdf= pdf + (ldt[1])*(ldt[2])*(1-ldt[3])*(ldt[4])*dmvnorm(x, mean = mean_x+mean_y1+mean_y2+mean_y4, sigma = cov_x+cov_y1+cov_y2+cov_y4)
  #1100
  pdf= pdf + (1-ldt[1])*(1-ldt[2])*(ldt[3])*(ldt[4])*dmvnorm(x, mean = mean_x+mean_y3+mean_y4, sigma = cov_x+cov_y3+cov_y4)
  #1101
  pdf= pdf + (ldt[1])*(1-ldt[2])*(ldt[3])*(ldt[4])*dmvnorm(x, mean = mean_x+mean_y1+mean_y3+mean_y4, sigma = cov_x+cov_y1+cov_y3+cov_y4)
  #1110
  pdf= pdf + (1-ldt[1])*(ldt[2])*(ldt[3])*(ldt[4])*dmvnorm(x, mean = mean_x+mean_y2+mean_y3+mean_y4, sigma = cov_x+cov_y2+cov_y3+cov_y4)
  #1111
  pdf= pdf + (ldt[1])*(ldt[2])*(ldt[3])*(ldt[4])*dmvnorm(x, mean = mean_x+mean_y1+mean_y2+mean_y3+mean_y4, sigma = cov_x+cov_y1+cov_y2+cov_y3+cov_y4)
  
  return(pdf)
}


negloglik_4assets_nocommon= function(params, x, dt, n) {
  # 
  # x is a matrix [Npoints * n] of all the points for which we compute the likelihood
  # 
  ## add check on inputs if necessary
  
  # reconstruction of parameters:
  idx =1
  mu = params[idx:(idx+n-1)]
  idx = idx+n 
  
  
  sigma=diag(params[idx:(idx+n-1)])
  idx = idx + n
  
  corr = matrix(rep(0,n*n), ncol = n)
  k=1
  for(i in 1:n)
    for(j in i:n){
      if (j!=i){
        corr[i,j]=params[idx+k-1]
        corr[j,i]=params[idx+k-1]
        k=k+1
      }
      else
        corr[i,j]=1
    }
  idx = idx + n*(n-1)/2
  
  S = sigma %*% corr %*% sigma
  
  # print(S)
  
  theta = params[idx:(idx+n-1)]
  idx = idx+n
  
  delta = params[idx:(idx+n-1)]
  idx = idx+n
  
  lambda = params[idx:(idx+n-1)]
  idx = idx+n
  
  # print(mu)
  # print(S)
  # print(theta)
  # print(delta)
  # print(lambda)
  # print(alpha)
  if( (idx-1)!=length(params))
    stop("Error in parameters reconstruction: number of parameters is wrong.")
  
  
  # computing pdf on each point and adding
  partial = MultivariateMertonPdf_4assets_nocommon(x, dt, mu, S, theta, delta, lambda)
  nll = -sum(log(partial))
  
  # last check on result
  if (is.nan(nll) | is.na(nll) | is.infinite(nll)) {
    nll = 1e10
  }
  return(nll)
}

