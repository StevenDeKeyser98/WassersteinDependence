# Penalized covariance matrix estimation using a penalty with derivative derpenal
# and nsteps iterations, default is 1
# Model selection is done by choosing penalty parameter in omega such that BIC is maximal
# S is ML estimator of the covariance matrix

covgpenal = function(S,n,omegas,derpenal = NULL,nsteps = 1){
  
  q = nrow(S) # Total dimension
  Sigma0 = S # Starting value
  delta = matrix(1,q,q) ; diag(delta) = 0 # Do not penalize diagonal 
  bic = integer(length(omegas)) # We want to maximize BIC
  ests = list() # Estimates
  
  # See equation (25) in paper
  
  for(o in 1:length(omegas)){
    
    omega = omegas[o]
    
    for(i in 1:nsteps){
      deltaT = delta * derpenal(delta * abs(Sigma0),omega) # Delta tilde matrix
      array = array(deltaT, dim = c(q,q,1))
      covgl = covglasso(S = S, n = n, lambda = array, start = S) # Solve (25)
      Sigma = covgl$sigma
      Sigma0 = Sigma # New starting value
    }
    bic[o] = as.numeric(covgl$BIC)
    ests[[o]] = covgl
  }
  best = which.max(bic)
  
  return(list(ests[[best]],omegas[best]))
}

SCAD = function(t,omega,a){
  
  # SCAD penalty
  
  if(t <= omega){
    return(omega*t)
    
  } else if(omega < t && t <= a*omega){
    
    return((1/(2*(a-1))) * (2*a*omega*t-t^2-omega^2))
  } else{
    
    return(((a+1)*(omega^2)/2))
    
  }
}

derSCAD = function(t,omega,a){
  
  # Derivative of SCAD penalty
  
  omega * ((t <= omega) + (max(a*omega - t,0)/((a-1) * omega)) * (t > omega))
  
}
