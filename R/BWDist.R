# Compute the Bures-Wasserstein distance between covariance matrices R and S

BWDist = function(R,S){ 
  
  sqrtR = Re(expm::sqrtm(R))
  BWD = Trace(R) + Trace(S) - 2*Trace(Re(expm::sqrtm(sqrtR %*% S %*% sqrtR)))
  
  return(BWD)
  
}
