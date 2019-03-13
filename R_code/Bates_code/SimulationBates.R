#################################################################################################
# Simulation of Multivariate Heston following 
# Dimitroff et al. "A Parsimonious Multi-Asset Heston Model: Calibration and Derivative Pricing"
#
# that we simply generalize to a SVJ Bates model by adding a jump component in the price dynamics
#
# NOTICE that we explicitly include the variance contribution to the drift in the log-returns
#################################################################################################
library(MASS)

simulate_mv_heston = function(Nasset=length(mu), mu, k,eta, theta, rho, CorrMatrix=matrix(1), 
                              S0=rep(1,Nasset), V0, dt=1/255, final_t=1, Nsim=100){
  
  if (Nasset!=length(mu) | Nasset!=length(k) | Nasset!=length(eta) | Nasset!=length(theta) 
      | Nasset!=length(rho) | Nasset!=length(S0) | Nasset!=length(V0) ){
    stop("Wrong dimension of parameters or initial values.")
  }
  
  full_corr = matrix(ncol = 2*Nasset, nrow = 2*Nasset)
  
  for (i in 1:Nasset){
    for (j in 1:Nasset) {
      full_corr[2*(i-1)+ c(1,2), 2*(j-1)+ c(1,2)] = correlation_block(rho,CorrMatrix, i,j)
    }
  }

  Nstep = ceil(final_t / dt)
  
  sim_x =array(0, dim = c(Nsim, Nstep+1, Nasset), dimnames = list(NULL,NULL,paste0("X_",1:Nasset)))
  sim_V = array(0, dim = c(Nsim, Nstep+1, Nasset),dimnames = list(NULL,NULL,paste0("V_",1:Nasset)))

  for (m in 1:Nasset) {
    sim_x[,1,m] = log(S0[m])
    sim_V[,1,m] = V0[m]
  }

  for (i in 1:(Nstep)) {
    z = mvrnorm(n = Nsim,mu=rep(0,Nasset*2),Sigma =full_corr)
    for (m in 1:Nasset) {
      
      V_plus= sim_V[,i,m]*(sim_V[,i,m]>0) # positive part for full truncation
      sim_V[,i+1, m] = sim_V[,i,m] + k[m]*(eta[m] - V_plus)*dt + theta[m]*sqrt(V_plus *dt)*z[,2*m]
      sim_x[, i+1,m] = sim_x[,i,m] + (mu[m]- V_plus*0.5)* dt + sqrt(V_plus * dt)*z[,2*m-1]
    }
  }
  
  res_V = lapply(seq(dim(sim_V)[3]), function(t) sim_V[ , , t])
  names(res_V) = paste0("V_",1:Nasset)
  res_x = lapply(seq(dim(sim_V)[3]), function(t) sim_x[ , , t])
  names(res_x) = paste0("x_",1:Nasset)
  return( list( return = res_x, variance = res_V ))
}




simulate_mv_bates = function(Nasset=length(mu), mu, k,eta, theta, rho,mu_j, sigma_j, lambda, CorrMatrix=matrix(1),
                             S0=rep(1,Nasset), V0, dt=1/255, final_t=1, Nsim=100){
  
  if (Nasset!=length(mu) | Nasset!=length(k) | Nasset!=length(eta) | Nasset!=length(theta) 
      | Nasset!=length(rho) | Nasset!=length(S0) | Nasset!=length(V0) |
      Nasset!=length(mu_j) | Nasset!=length(sigma_j) | Nasset!=length(lambda)){
    stop("Wrong dimension of parameters or initial values.")
  }
  
  # correlation for the brownian motions driving (S1,V1,S2,V2, ... ,Sn,Vn)
  full_corr = matrix(ncol = 2*Nasset, nrow = 2*Nasset) 
  for (i in 1:Nasset){
    for (j in 1:Nasset) {
      full_corr[2*(i-1)+ c(1,2), 2*(j-1)+ c(1,2)] = correlation_block(rho,CorrMatrix, i,j)
    }
  }
  
  full_corr=regularization(full_corr)
  
  Nstep = ceil(final_t / dt)
  
  sim_x =array(0, dim = c(Nsim, Nstep+1, Nasset), dimnames = list(NULL,NULL,paste0("X_",1:Nasset)))
  sim_V = array(0, dim = c(Nsim, Nstep+1, Nasset),dimnames = list(NULL,NULL,paste0("V_",1:Nasset)))
  
  for (m in 1:Nasset) {
    sim_x[,1,m] = log(S0[m])
    sim_V[,1,m] = V0[m]
  }
  
  for (i in 1:(Nstep)) {
    z = mvrnorm(n = Nsim,mu=rep(0,Nasset*2) ,Sigma =full_corr)
    for (m in 1:Nasset) {
      dq = rbinom(Nsim,size = 1, prob= lambda[m]*dt)
      jump = rnorm(Nsim, mean = log(1+mu_j[m]) - 0.5*sigma_j[m]^2, sd = sigma_j[m])
      V_plus= sim_V[,i,m]*(sim_V[,i,m]>0) # positive part for full truncation
      sim_V[,i+1, m] = sim_V[,i,m] + k[m]*(eta[m] - V_plus)*dt + theta[m]*sqrt(V_plus *dt)*z[,2*m]
      sim_x[, i+1,m] = sim_x[,i,m] + (mu[m] - V_plus*0.5 -lambda[m]*mu_j[m] )* dt + sqrt(V_plus * dt)*z[,2*m-1] + dq * (jump)
    }
  }
  
  res_V = lapply(seq(dim(sim_V)[3]), function(t) sim_V[ , , t])
  names(res_V) = paste0("V_",1:Nasset)
  res_x = lapply(seq(dim(sim_V)[3]), function(t) sim_x[ , , t])
  names(res_x) = paste0("x_",1:Nasset)
  return( list( return = res_x, variance = res_V ))
}



# Aux function to compute correlation block as in the reference paper
correlation_block= function(rho, CorrMatrix, i, j){
  if (i==j){
    res = rbind(c(1,  rho[i]),
                c(rho[i],  1))
  }
  else{
    res = CorrMatrix[i,j] * rbind( c(1, rho[j]), c(rho[i], rho[i]*rho[j])) 
  }
  res
}




# computes correlation matrix from given simulation result
correlation_MC_estimation = function( simulation_result){
  # simulation result should be a n-asset simulation output of simulate_mv_heston or bates
  # NOTE: correlation is computed on the log returns in [t,t+dt] 
  
  Nasset = length(simulation_result$return) 
  Nsim = dim(simulation_result$return[[1]])[1] 
  
  if( Nasset == 2 ){
    sumcorr= 0
    for (k in 1:Nsim) {
      # correlation of log-returns in path k
      sumcorr = sumcorr + cor(diff(simulation_result$return[[1]][k,]), diff(simulation_result$return[[2]][k,]))
    }
    res = sumcorr/Nsim      # average correlation
  }
  else if (Nasset > 2){
    corr_matrix = diag(nrow = Nasset) # Identity matrix 
    for ( i in 1:(Nasset-1)) {
      for (j in (i+1):Nasset) {
        sumcorr = 0
        for (k in 1:Nsim) {
          sumcorr = sumcorr + cor(diff(simulation_result$return[[i]][k,]), diff(simulation_result$return[[j]][k,]))
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




# Computes correlation matrix from initial and model parameters [ used in model correlation calibration]
expected_model_correlation = function(model_corr, mu, k,eta, theta, rho,mu_j, sigma_j, lambda, S0, V0, dt =1/255, final_t = 1, Nsim){
  simulated = simulate_mv_bates(Nasset = 2, mu, k,eta, theta, rho, mu_j, sigma_j, lambda, matrix(c(1,model_corr, model_corr,1),ncol=2),
                                S0,V0, dt, final_t, Nsim = Nsim)
  rho_MC = correlation_MC_estimation(simulated)
  rho_MC
}




# calibrates the model correlation of a 2-asset heston or bates given the sample corr
calibrate_correlation = function(sample_corr,Nasset=length(mu), mu, k,eta, theta, rho,mu_j, sigma_j, lambda, S0=c(1,1), V0,  dt=1/255, final_t=1, Nsim=100){
  if (Nasset!=length(mu) | Nasset!=length(k) | Nasset!=length(eta) | Nasset!=length(theta) 
      | Nasset!=length(rho) | Nasset!=length(S0) | Nasset!=length(V0) |
      Nasset!=length(mu_j) | Nasset!=length(sigma_j) | Nasset!=length(lambda)){
    stop("Wrong dimension of parameters or initial values.")
  }
  
  f_to_be_solved = function(model_corr,  sample_corr, ...){
    expected_model_correlation(model_corr=model_corr, ...) - sample_corr
  }
  
  # print(Nsim)
  # print(test)
  res = uniroot(f=f_to_be_solved, interval =c(-1,1),
                sample_corr = sample_corr,
                mu = mu, k=k, eta=eta, theta=theta, rho=rho, mu_j=mu_j, sigma_j=sigma_j, lambda = lambda,
                S0=S0, V0=V0, dt=dt, final_t=final_t, Nsim=Nsim,
                trace = 2)
  
  return(res$root)
}


# Calibrates the n*n correlation matrix for a n-asset model by iterating on all pairs of assets
calibrate_full_correlation = function(sample_corr_matrix,Nasset=length(mu), mu, k,eta, theta, rho,mu_j=rep(0,Nasset), sigma_j=rep(1,Nasset), lambda=rep(0,Nasset), 
                                      S0=rep(1,Nasset), V0,  dt=1/255, final_t=1, Nsim=100){
  
  if (Nasset!=length(mu) | Nasset!=length(k) | Nasset!=length(eta) | Nasset!=length(theta) 
      | Nasset!=length(rho) | Nasset!=length(S0) | Nasset!=length(V0) |
      Nasset!=length(mu_j) | Nasset!=length(sigma_j) | Nasset!=length(lambda)){
    stop("Wrong dimension of parameters or initial values.")
  }
  
  model_corr_matrix = diag(Nasset)
  #print(model_corr_matrix)
  
  for (i in 1:(Nasset-1)) {
    for (j in (i+1):Nasset) {
      print(paste0("Estimating element [", i, ',', j,']'))
      
      upper_val = expected_model_correlation(model_corr = 1, mu=mu[c(i,j)], k=k[c(i,j)],eta=eta[c(i,j)], 
                                             theta=theta[c(i,j)], rho=rho[c(i,j)], 
                                             mu_j=mu_j[c(i,j)], sigma_j=sigma_j[c(i,j)], lambda=lambda[c(i,j)], S0=c(1,1),
                                             V0=V0[c(i,j)], dt=dt, final_t=final_t, Nsim=Nsim) -sample_corr_matrix[i,j]
      
      lower_val = expected_model_correlation(model_corr = -1, mu=mu[c(i,j)], k=k[c(i,j)],eta=eta[c(i,j)], 
                                             theta=theta[c(i,j)], rho=rho[c(i,j)], 
                                             mu_j=mu_j[c(i,j)], sigma_j=sigma_j[c(i,j)], lambda=lambda[c(i,j)],  S0=c(1,1),
                                             V0=V0[c(i,j)], dt=dt, final_t=final_t, Nsim=Nsim) -sample_corr_matrix[i,j]
      
      if(upper_val*lower_val < 0){
        model_corr_matrix[i,j] =  calibrate_correlation(sample_corr_matrix[i,j], mu=mu[c(i,j)], k=k[c(i,j)],eta=eta[c(i,j)], 
                                                        theta=theta[c(i,j)], rho=rho[c(i,j)], 
                                                        mu_j=mu_j[c(i,j)], sigma_j=sigma_j[c(i,j)], lambda=lambda[c(i,j)],
                                                        V0=V0[c(i,j)], dt=dt, final_t=final_t, Nsim=Nsim)
      }
      else{
        print(paste("Cannot obtain a correlation of", sample_corr_matrix[i,j] , "with given parameters." ))
        if(abs(upper_val)>abs(lower_val)){
          model_corr_matrix[i,j] = -10
        }
        else{
          model_corr_matrix[i,j] = 10
        }
      }
      
      
      model_corr_matrix[j,i] = model_corr_matrix[i,j] 
      # print(model_corr_matrix)
    }
  }
  
  # print(model_corr_matrix)
  #model_corr_matrix = regularization(model_corr_matrix)
  
  return(model_corr_matrix)
}


# Function to regularize the mode correlation matrix in case the resulting one 
# is not positive definite [in the end we used a trivial procedure as default rather than jackel
# because the latter is less reliable and might produce slightly negative eigenvalues e.g. -1e-6]
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





