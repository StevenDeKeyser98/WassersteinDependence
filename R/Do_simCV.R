  # Simulation for covglasso on correlation matrix R using samples of sizes in vector n
  # and M Monte Carlo replications. Either use the adaptive lasso with rho values
  # or equal penalization to all off diagonal elements with omega values.
  # Group lasso is also possible.
  # Best model gets selected by 5-fold cross-validation.



Do_simCV = function(R,dim,n,M,rho = NULL, omega = NULL, alpha = NULL, scad = FALSE, group = FALSE){
  
  q = nrow(R)
  R_ests = list() # ML estimates 
  
  if(!is.null(omega)){
    est_omega = list()  # Selected omegas
    P = matrix(1,q,q) # Penalize all off-diagonal entries
    diag(P) = 0 # Do not penalize diagonal elements
  } 
  
  if(!is.null(rho)){est_rho = list()} # Selected rhos 
  
  for(i in 1:length(n)){
    R_ests[[paste("n =", toString(n[i]))]] = list()
    if(!is.null(omega)){est_omega[[paste("n =", toString(n[i]))]] = integer(M)}
    if(!is.null(rho)){est_rho[[paste("n =", toString(n[i]))]] = integer(M)}
  }
  
  R_ests_covglasso = R_ests # Covglasso estimates 
  
  set.seed(123)
  for(i in 1:length(n)){
    for(e in 1:M){
      print(paste("n = ",n[i]," : Monte Carlo sample",e,"out of",M))
      Gsample = rmvnorm(n[i],rep(0,q),R, method = "chol") # Multivariate Gaussian sample
      Usample = pnorm(Gsample) # Multivariate Gaussian copula sample, marginals can be changed
      sample = qnorm(Usample) # We take standard normal marginals for now
      
      scores = matrix(0,n[i],q)
      for(j in 1:q){
        scores[,j] = qnorm((n[i]/(n[i]+1)) * ecdf(sample[,j])(sample[,j]))
      }
      
      Sigma_est = cov(scores) * ((n[i]-1)/n[i])
      
      R_est = est_R(sample,1)
      
      # We now estimate the covariance matrix of the pseudo Gaussian sample using a Lasso penalty (package covglasso)
      
      if(!is.null(omega) && scad == FALSE && group == FALSE){
        est_rho = NULL
        best_omega = CVomega2(scores, k = 5, omegas = omega, method = "Lasso")
        covglasso = covglasso(S = Sigma_est, n = n[i], lambda = array(best_omega*P, dim = c(q,q,1)), start = Sigma_est)
        est_omega[[paste("n =", toString(n[i]))]][e] = best_omega
      }
      
      if(!is.null(rho)){
        est_omega = NULL
        best_rho = CVomega2(scores, k = 5, rhos = rho, method = "ALasso")
        covglasso = covglasso(S = Sigma_est, n = n[i],  rho = best_rho, start = Sigma_est)
        est_rho[[paste("n =", toString(n[i]))]][e] = best_rho
      }
      
      if(scad == TRUE){
        est_rho = NULL
        best_omega  = CVomega2(scores, k = 5, omegas = omega, method = "Scad")
        PT = P * apply(P * abs(Sigma_est),1:2,function(t){derSCAD(t,best_omega,3.7)})
        covglasso = covglasso(S = Sigma_est, n = n[i], lambda = array(PT, dim = c(q,q,1)), start = Sigma_est)
        est_omega[[paste("n =", toString(n[i]))]][e] = best_omega
      }
      
      if(group == TRUE){
        est_rho = NULL
        best_omega = CVomega2(scores, k = 5, dim = dim, omegas = omega, method = "GLasso")
        covglasso = Gspcov(Sigma_est, Sigma_est, best_omega, dim, step.size = 100)
        names(covglasso)[[2]] = "sigma"
        est_omega[[paste("n =", toString(n[i]))]][e] = best_omega
      }
      
      Sigma_est_covglasso = covglasso$sigma
      D = sqrt(solve(diag(diag(Sigma_est_covglasso))))
      R_est_covglasso = D %*% Sigma_est_covglasso %*% D
      
      R_ests[[paste("n =", toString(n[i]))]][[e]] = R_est
      R_ests_covglasso[[paste("n =", toString(n[i]))]][[e]] = R_est_covglasso
    }
  }
  
  FerrorsR_est = list() # Frobenius errors ML estimates
  for(i in 1:length(n)){
    FerrorsR_est[[paste("n =", toString(n[i]))]] = integer(M)
  }
  
  FerrorsR_est_covglasso = FerrorsR_est # Frobenius errors covglasso estimates
  TPR = FerrorsR_est # True positive rates
  FPR = FerrorsR_est # False positive rates
  MIests = FerrorsR_est # Mutual information estimates for ML estimates 
  Helests = FerrorsR_est # Hellinger distance estimates for ML estimates 
  BWD1ests = FerrorsR_est # BWD1 estimates for ML estimates
  BWD2ests = FerrorsR_est # BWD2 estimates for ML estimates 
  MIests_covglasso = FerrorsR_est # Mutual information estimates for covglasso estimates 
  Helests_covglasso = FerrorsR_est # Hellinger distance estimates for covglasso estimates
  BWD1ests_covglasso = FerrorsR_est # BWD1 estimates for covglasso estimates
  BWD2ests_covglasso = FerrorsR_est # BWD2 estimates for covglasso estimates
  
  
  for(i in 1:length(n)){
    for(e in 1:M){
      R_est = R_ests[[paste("n =", toString(n[i]))]][[e]]
      R_est_covglasso = R_ests_covglasso[[paste("n =", toString(n[i]))]][[e]]
      TP = length(which(as.numeric(R_est_covglasso == 0) + as.numeric(R == 0) == 2)) # True positive
      AP = sum(as.numeric(R == 0)) # Actual positives
      FP = length(which(as.numeric(R_est_covglasso == 0) + as.numeric(R != 0) == 2)) # False positive
      TPR[[paste("n =", toString(n[i]))]][e] = TP/AP
      FPR[[paste("n =", toString(n[i]))]][e] = FP/((q^2)-AP)
      FerrorsR_est[[paste("n =", toString(n[i]))]][e] = sqrt(mean((R - R_est)^2))
      FerrorsR_est_covglasso[[paste("n =", toString(n[i]))]][e] = sqrt(mean((R - R_est_covglasso)^2))
      MIests[[paste("n =", toString(n[i]))]][e] = MI_normal(R_est,dim)
      Helests[[paste("n =", toString(n[i]))]][e] = Hel_normal(R_est,dim)
      BWD1ests[[paste("n =", toString(n[i]))]][e] = BWD1(R_est,dim)
      BWD2ests[[paste("n =", toString(n[i]))]][e] = BWD2(R_est,dim)
      MIests_covglasso[[paste("n =", toString(n[i]))]][e] = MI_normal(R_est_covglasso,dim)
      Helests_covglasso[[paste("n =", toString(n[i]))]][e] = Hel_normal(R_est_covglasso,dim)
      BWD1ests_covglasso[[paste("n =", toString(n[i]))]][e] = BWD1(R_est_covglasso,dim)
      BWD2ests_covglasso[[paste("n =", toString(n[i]))]][e] = BWD2(R_est_covglasso,dim)
    }
  }
  
  return(list("R_ests" = R_ests, "R_ests_covglasso" = R_ests_covglasso,
              "TPR" = TPR, "FPR" = FPR, "est_omega" = est_omega, "est_rho" = est_rho,
              "FerrorsR_est" = FerrorsR_est, "FerrorsR_est_covglasso" = FerrorsR_est_covglasso, 
              "MIests" = MIests, "Helests" = Helests, "BWD1ests" = BWD1ests, "BWD2ests" = BWD2ests,
              "MIests_covglasso" = MIests_covglasso, "Helests_covglasso" = Helests_covglasso,
              "BWD1ests_covglasso" = BWD1ests_covglasso, "BWD2ests_covglasso" = BWD2ests_covglasso))
}

# Function for cross-validation

CVomega2 = function(data,k,dim = NULL,omegas = NULL,rhos = NULL,method){
  
  # Selects omega or rho for sparse covariance matrix estimation via k fold cross validation
  # Data should be the matrix of sample normal scores 
  # Omegas is the grid of possible omega values
  # rhos is the grid of possible rho values
  # method is the penalization technique used
  
  q = ncol(data)
  P = matrix(1,q,q) # Penalize all off-diagonal entries
  diag(P) = 0 # Do not penalize diagonal elements
  n = nrow(data) ; flds = list() ; l = floor(n/k)
  for(i in 1:k){flds[[i]] = ((i-1)*l + 1) : (i*l)}
  if(n > (k*l)){for(i in 1:(n-k*l)){flds[[i]] = c(flds[[i]],k*l+i)}}
  
  if(!is.null(omegas)){
    CVLL = integer(length(omegas))
    for(i in 1:length(omegas)){
      omega = omegas[i]
      Cvomega = integer(k)
      for(j in 1:k){
        test = data[flds[[j]],]
        train = data[-flds[[j]],]
        cov_test = cov(test)
        cov_train = cov(train)
        nj = nrow(train)
        if(method == "Lasso"){
          cov_train_penal = covglasso(S = cov_train, n = nj, lambda = array(omega*P, dim = c(q,q,1)), start = cov_train)$sigma
        }
        if(method == "Scad"){
          PT = P * apply(P * abs(cov_train),1:2,function(t){derSCAD(t,omega,3.7)})
          cov_train_penal = covglasso(S = cov_train, n = nj, lambda = array(PT, dim = c(q,q,1)), start = cov_train)$sigma
        }
        if(method == "GLasso"){
          cov_train_penal = Gspcov(cov_train, cov_train, omega, dim, step.size = 100)$Sigma
        }
        Cvomega[j] = log(det(cov_train_penal)) + sum(diag(cov_test %*% solve(cov_train_penal)))
      }
      CVLL[i] = mean(Cvomega)
    }
    opt_omega = omegas[which.min(CVLL)]
    return(opt_omega)
  }
  
  if(!is.null(rhos)){
    CVLL = integer(length(rhos))
    for(i in 1:length(rhos)){
      rho = rhos[i]
      Cvrho = integer(k)
      for(j in 1:k){
        test = data[flds[[j]],]
        train = data[-flds[[j]],]
        cov_test = cov(test)
        cov_train = cov(train)
        nj = nrow(train)
        cov_train_penal = covglasso(S = cov_train, n = nj,  rho = rho, start = cov_train)$sigma
        Cvrho[j] = log(det(cov_train_penal)) + sum(diag(cov_test %*% solve(cov_train_penal)))
      }
      CVLL[i] = mean(Cvrho)
    }
    
    opt_rho = rhos[which.min(CVLL)]
    return(opt_rho)
    
  }
}
