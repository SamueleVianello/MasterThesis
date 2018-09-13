ProcessoMerton_Correlazione <-function(mat_cholesky, drift,sigma,N,lambdaj,mu_xi,sigma_q,h,n_asset){
  
  #restituisce una matrice con l'evoluzione di 3 asset con matrice di correlazione assegnata
  # e anche i salti nel processo di Poisson contengono una correlazione
  
  
  
  t=(0:N/h)/N
  X=matrix(nrow=n_asset, ncol= N+1)
  F=matrix(nrow=n_asset, ncol= N+1)
  I=matrix(nrow=n_asset,ncol=N) # I contiene i salti di Asset1, Asset2, Asset3, poi i salti in comune tra 1 e 2, 1 e 3, 2 e 3,
  # la settima riga ha i crolli in comune tra tutti gli asset
  
  
  X[,1]=c(1,1,1)  # inizializzo a 1 tutti gli asset
  
  for(j in 1:N) {  #ciclo sui mesi, da 1 a 120
    # consideriamo 70% di salti individuali, 20% di correlazione a due a due, 10% di salti in comune

    individuale1=PPgen(h*lambdaj[1])  
    individuale2=PPgen(h*lambdaj[2])  
    individuale3=PPgen(h*lambdaj[3])  
    I[1,j]=individuale1
    I[2,j]=individuale2
    I[3,j]=individuale3
    
    # #browser()
    # 
    # comune_12=PPgen(h*lambdaj*0.1)        #   e' stata tolta la correlazione sui salti, reintrodurla se necessario
    # comune_13=PPgen(h*lambdaj*0.1) 
    # comune_23=PPgen(h*lambdaj*0.1) 
    # 
    # comune_123=PPgen(h*lambdaj*0.1)
    # 
    # I[1,j]=individuale1+ comune_12+ comune_13+ comune_123
    # I[2,j]=individuale2+ comune_12+ comune_23+ comune_123
    # I[3,j]=individuale3+ comune_13+ comune_23+ comune_123
    # # browser()
    
    if (I[1,j]==0){F[1,j]=0} else {
      F[1,j]=rnorm(1, mean=mu_xi[1],sd=sigma_q[1]) 
    }
    if (I[2,j]==0){F[2,j]=0} else {
      F[2,j]=rnorm(1, mean=mu_xi[2],sd=sigma_q[2])  
    }
    if (I[3,j]==0){F[3,j]=0} else {
      F[3,j]=rnorm(1, mean=mu_xi[3],sd=sigma_q[3])  
    }
    
    gaussiana_correlazione = mat_cholesky %*% rnorm(n_asset)
    
    
    # X[,j+1]=X[,j] + drift*h*X[,j]+sigma*sqrt(h)*gaussiana_correlazione*X[,j]-F[,j]*X[,j]
    X[,j+1]=X[,j]*(1 + drift*h +sigma*sqrt(h)*gaussiana_correlazione-F[,j])
    
  }
  # browser()
  return (X)
}