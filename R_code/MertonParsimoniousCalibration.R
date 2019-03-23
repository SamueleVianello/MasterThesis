#parsimonious mvmerton

library(MASS)

simulate_mv_merton = function(Nasset=length(mu), mu, vol=NA, mu_j, sigma_j, lambda, CorrMatrix=matrix(NA), CovMatrix = matrix(NA),
                              x0=rep(0,Nasset), dt=1/255, final_t=1, Nsim=100){
  
  if(is.na(CorrMatrix[1])& Nasset==1){
    CorrMatrix=matrix(1)
  }
  if(is.na(CovMatrix[1])){
    if (is.na(vol[1]) | is.na(CorrMatrix[1])){
      stop("Need to specify volatilities and Correlation Matrix, or only Covariance Matrix")
    }
  }
  if(is.na(vol[1]) & is.na(CorrMatrix[1])){
    CorrMatrix = cov2cor(CovMatrix)
    vol = sqrt(diag(CovMatrix))
  }
  
  
  if (Nasset!=length(mu) | Nasset!=length(x0) |
      Nasset!=length(mu_j) | Nasset!=length(sigma_j) | Nasset!=length(lambda)){
    #print(paste(Nasset,length(mu),length(x0),length(mu_j),length(sigma_j),length(lambda)))
    stop("Wrong dimension of parameters or initial values.")
  }
  
  Nstep = ceil(final_t / dt)
  sim_x =array(0, dim = c(Nsim, Nstep+1, Nasset), dimnames = list(NULL,NULL,paste0("X_",1:Nasset)))
  for (m in 1:Nasset) {
    sim_x[,1,m] = x0[m]
  }

  for (i in 1:(Nstep)) {
    z = mvrnorm(n = Nsim,mu=rep(0,Nasset) ,Sigma =CorrMatrix)
    for (m in 1:Nasset) {
      dq = rpois(Nsim, lambda = lambda[m]*dt)
      # jump = rnorm(Nsim, mean = log(1+mu_j[m]) - 0.5*sigma_j[m]^2, sd = sigma_j[m])
      jump = rnorm(Nsim, mean = mu_j[m], sd = sigma_j[m])
      sim_x[, i+1,m] = sim_x[,i,m] + (mu[m] - vol[m]^2*0.5 -lambda[m]*mu_j[m] )* dt + vol[m]*sqrt(dt)*z[,m] + dq * (jump)
    }
  }
  
  res_x = lapply(seq(dim(sim_x)[3]), function(t) sim_x[ , , t])
  names(res_x) = paste0("x_",1:Nasset)
  return( res_x)
}



# computes correlation matrix from given simulation result
correlation_MC_estimation = function( simulation_result){
  # simulation result should be a n-asset simulation output of simulate_mv_merton
  # NOTE: correlation is computed on the log returns in [t,t+dt] 
  
  Nasset = length(simulation_result) 
  Nsim = dim(simulation_result[[1]])[1] 
  
  if( Nasset == 2 ){
    sumcorr= 0
    for (k in 1:Nsim) {
      # correlation of log-returns in path k
      sumcorr = sumcorr + cor(diff(simulation_result[[1]][k,]), diff(simulation_result[[2]][k,]))
    }
    res = sumcorr/Nsim      # average correlation on different scenarios
  }
  else if (Nasset > 2){
    corr_matrix = diag(nrow = Nasset) # Identity matrix 
    for ( i in 1:(Nasset-1)) {
      for (j in (i+1):Nasset) {
        sumcorr = 0
        for (k in 1:Nsim) {
          sumcorr = sumcorr + cor(diff(simulation_result[[i]][k,]), diff(simulation_result[[j]][k,]))
        }
        corr_matrix[i,j] = sumcorr/Nsim
        corr_matrix[j,i] = sumcorr/Nsim
      }
    }
    res = corr_matrix
  }
  else{
    stop("Need at least 2 assets to compute correlation.")
  }
  return(res)
}



# Computes correlation for 2 assets from initial and model parameters [ used in model correlation calibration]
expected_model_correlation_merton = function(model_corr, mu, vol, mu_j, sigma_j, lambda,
                                      x0, dt=1/255, final_t=1, Nsim){
  simulated = simulate_mv_merton(Nasset = 2, mu, vol, mu_j, sigma_j, lambda, matrix(c(1,model_corr, model_corr,1),ncol=2),
                                x0=x0, dt=dt, final_t=final_t, Nsim = Nsim)
  rho_MC = correlation_MC_estimation(simulated)
  rho_MC
}


# calibrates the model correlation of a 2-asset heston or bates given the sample corr
calibrate_correlation_merton = function(sample_corr,Nasset=length(mu), mu, vol, mu_j, sigma_j, lambda, x0=c(0,0), dt=1/255, final_t=1, Nsim=100){
  
  if (Nasset!=length(mu) | Nasset!=length(x0) |
      Nasset!=length(mu_j) | Nasset!=length(sigma_j) | Nasset!=length(lambda)){
    #print(paste(Nasset,length(mu),length(x0),length(mu_j),length(sigma_j),length(lambda)))
    stop("Wrong dimension of parameters or initial values.")
  }
  
  f_to_be_solved = function(model_corr,  sample_corr, ...){
    expected_model_correlation_merton(model_corr=model_corr, ...) - sample_corr
  }
  
  # print(Nsim)
  # print(test)
  res = uniroot(f=f_to_be_solved, interval =c(-1,1),
                sample_corr = sample_corr,
                mu = mu, vol=vol, mu_j=mu_j, sigma_j=sigma_j, lambda = lambda,
                x0=x0, dt=dt, final_t=final_t, Nsim=Nsim,
                trace = 3)
  
  return(res$root)
}



# Calibrates the n*n correlation matrix for a n-asset model by iterating on all pairs of assets
calibrate_full_correlation_merton = function(sample_corr_matrix, Nasset=length(mu), mu, vol, mu_j, sigma_j, lambda, 
                                      x0=rep(0,Nasset), dt=1/255, final_t=1, Nsim=500){
  
  if (Nasset!=length(mu) | Nasset!=length(x0) |
      Nasset!=length(mu_j) | Nasset!=length(sigma_j) | Nasset!=length(lambda)){
    #print(paste(Nasset,length(mu),length(x0),length(mu_j),length(sigma_j),length(lambda)))
    stop("Wrong dimension of parameters or initial values.")
  }
  
  # create identity matrix
  model_corr_matrix = diag(Nasset) 
  
  for (i in 1:(Nasset-1)) {
    for (j in (i+1):Nasset) {
      print(paste0("Estimating element [", i, ',', j,']'))
      upper_val = expected_model_correlation_merton(model_corr = 1, mu=mu[c(i,j)], vol = vol[c(i,j)],
                                                    mu_j=mu_j[c(i,j)], sigma_j=sigma_j[c(i,j)], lambda=lambda[c(i,j)], 
                                                    x0=x0[c(i,j)], dt=dt, final_t=final_t, Nsim=Nsim) -sample_corr_matrix[i,j]
      lower_val = expected_model_correlation_merton(model_corr = -1, mu=mu[c(i,j)], vol = vol[c(i,j)],
                                                    mu_j=mu_j[c(i,j)], sigma_j=sigma_j[c(i,j)], lambda=lambda[c(i,j)], 
                                                    x0=x0[c(i,j)], dt=dt, final_t=final_t, Nsim=Nsim) -sample_corr_matrix[i,j]
      if(upper_val*lower_val < 0){
      model_corr_matrix[i,j] =  calibrate_correlation_merton(sample_corr_matrix[i,j], mu=mu[c(i,j)], vol = vol[c(i,j)],
                                                      mu_j=mu_j[c(i,j)], sigma_j=sigma_j[c(i,j)], lambda=lambda[c(i,j)], 
                                                      x0=x0[c(i,j)], dt=dt, final_t=final_t, Nsim=Nsim)
      }
      else{
        print(paste("Cannot obtain a correlation of", sample_corr_matrix[i,j] , "with given parameters." ))
        if(abs(upper_val)>abs(lower_val)){
          model_corr_matrix[i,j] = -1
        }
        else{
          model_corr_matrix[i,j] = 1
        }
      }
      model_corr_matrix[j,i] = model_corr_matrix[i,j]
      # print(model_corr_matrix)
    }
    #print(model_corr_matrix[i,])
  }
  
  # print(model_corr_matrix)
  
 model_corr_matrix = regularization(model_corr_matrix, method = "jackel")
  
  return(model_corr_matrix)
}

regularization = function(mat, method =  'simple'){
  decomp = eigen(mat,symmetric = TRUE)
  S = decomp$vectors
  eigval = decomp$values
  
  # print("Eigenvalues: ")
  # print(eigval)
  
  if(method =="jackel"){
    eigval[which(eigval<0)]=0
    
    Lambda = diag(eigval)
    
    t = rep(x=NA, length(eigval))
    for (i in 1:length(eigval)) {
      t[i] = 1/sum(S[i,]^2 * eigval)
    }  
    B = sqrt(diag(t)) %*% S %*% sqrt(Lambda)
    
    res = B %*% t(B)
  }
  else{
    eigval[which(eigval<0)]=1e-5
    res = S %*% diag(eigval)%*% t(S)
  }
  res
}

