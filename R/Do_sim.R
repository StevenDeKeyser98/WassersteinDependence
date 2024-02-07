# Simulation for covglasso on correlation matrix R using samples of sizes in vector n
# and M Monte Carlo replications. Either use the adaptive lasso with rho values
# or equal penalization to all off diagonal elements with omega values 
# Best model gets selected by BIC 
# Dim contains the dimensions of the random vectors

Do_sim = function(R,dim,n,M,rho = NULL, omega = NULL, scad = FALSE, group = FALSE){
  
  q = nrow(R)
  R_ests = list() # ML estimates 
  
  if(!is.null(omega)){
    est_omega = list()  # Selected omegas
    P = matrix(1,q,q) # Penalize all off-diagonal entries
    diag(P) = 0 # Do not penalize diagonal elements
    array = array(0,dim = c(q,q,length(omega)))
    
    for(o in 1:length(omega)){array[,,o] = omega[o] * P}
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
      sample = qnorm(Usample) # We take standard normal marginals 
      
      scores = matrix(0,n[i],q) # Normal scores
      
      for(j in 1:q){
        scores[,j] = qnorm((n[i]/(n[i]+1)) * ecdf(sample[,j])(sample[,j]))
      }
      
      Sigma_est = cov(scores) * ((n[i]-1)/n[i]) # ML estimator for sigma
      
      R_est = est_R(sample,1)
      
      # We now estimate the covariance matrix of the pseudo Gaussian sample using a Lasso penalty (package covglasso)
      
      if(!is.null(omega) && scad == FALSE && group == FALSE){ # Lasso
        est_rho = NULL
        covglasso = covglasso(S = Sigma_est, n = n[i], lambda = array, start = Sigma_est)
        est_omega[[paste("n =", toString(n[i]))]][e] = omega[which.max(covglasso$BIC)]
      }
      
      if(!is.null(rho)){ # Adaptive lasso
        est_omega = NULL
        covglasso = covglasso(S = Sigma_est, n = n[i],  rho = rho, start = Sigma_est)
        est_rho[[paste("n =", toString(n[i]))]][e] = rho[which.max(covglasso$BIC)]
      }
      
      if(scad == TRUE){ # Scad 
        est_rho = NULL
        covgpen = covgpenal(Sigma_est,n[i],omega,derpenal = function(t,omega){derSCAD(t,omega,3.7)})
        covglasso = covgpen[[1]]
        est_omega[[paste("n =", toString(n[i]))]][e] = covgpen[[2]]
      }
      
      if(group == TRUE){ # Group lasso
        est_rho = NULL
        groupl = GroupLasso(Sigma_est, Sigma_est, n[i], omega, dim, step.size = 100)
        covglasso = groupl[[1]]
        est_omega[[paste("n =", toString(n[i]))]][e] = groupl[[2]]
      }
      
      Sigma_est_covglasso = covglasso$sigma # Penalized sigma estimate
      D = sqrt(solve(diag(diag(Sigma_est_covglasso)))) 
      R_est_covglasso = D %*% Sigma_est_covglasso %*% D # Penalized R estimate
      
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
