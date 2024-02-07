# Find penalty parameter that maximizes Gaussian cross-validated 
# likelihood using K folds
# Candidate values are contained in omegas 
# Data is a matrix with the observations in the rows

CVomega = function(data,omegas,K){
  
  n = nrow(data) # Sample size
  s = ncol(data) # Total dimension
  scores = matrix(0,n,s) # Normal scores
  
  for(j in 1:s){ # Compute normal scores
    scores[,j] = qnorm((n/(n+1)) * ecdf(data[,j])(data[,j]))
  }
  
  LLs = integer(length(omegas))
  
  for(l in 1:length(omegas)){ # Compute CVLF for every candidate penatly parameter
    LLs[l] = CVLF(omegas[l],scores,K)
  }
  
  return(omegas[which.max(LLs)])
  
}


LogL = function(data,R,omega){
  
  # Log-likelihood of Gaussian model with correlation matrix R and penalty parameter omega
  # Data is a matrix with the observations in the rows
  
  s = nrow(R) # Total dimension
  n = nrow(data) # Sample size
  
  # See equation (5) in Warton (2008)
  
  sigmaD = diag(rep(sqrt(sum(qnorm(seq(1,n)/(n+1))^2)/(n-1)),s))
  sigmaL = (1/omega) * sigmaD %*% (omega * R + (1-omega) * diag(s)) %*% sigmaD
  t1 = n * s * log(2*pi) + n * log(det(sigmaL))
  t2 = sum(diag(data %*% solve(sigmaL) %*% t(data)))
  
  return((-1/2)*(t1+t2))
  
}

CVLF = function(omega,data,K){
  
  # Cross-validated likelihood function with penalty parameter omega,
  # and K folds for Ridge penalization
  # Data is a matrix with the normal scores in the rows
  
  n = nrow(data) # Sample size
  ind = sample(seq(1,n)) # Random assignment 
  fld_size = floor(n/K) # Size of each fold
  flds = list() # Folds
  
  for(i in 1:K){ # Construct folds
    flds[[i]] = ind[((fld_size)*(i-1) + 1):(fld_size*i)]
  }
  
  if(n-(K*fld_size) > 0){
    # In case n is not a multiple of K, assign remaining elements to first folds
    for(i in 1:(n-(K*fld_size))){
      flds[[i]] = c(flds[[i]],ind[fld_size*K + i])
    }
  }
  LL = 0 # Cross-validated likelihood, see equation (4) in Warton (2008)
  
  for(i in 1:K){
    valid = data[flds[[i]],] # Validation data
    train = data[setdiff(seq(1,n),flds[[i]]),] # Training data
    R_est = cor(train)
    LL = LL + LogL(valid,R_est,omega)
  }
  
  return(LL)
  
}
