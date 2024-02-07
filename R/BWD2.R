# Computes the second BW dependence coefficient given the correlation matrix R
# of the random vectors of dimensions dim

BWD2 = function(R,dim){
  
  s = nrow(R) # Total dimension
  k = length(dim) # Amount of random vectors 
  R0 = create_R0(R,dim) # R0
  eigen_Rii = matrix(0,s,k) # Eigenvalues of diagonal blocks
  start = 1 # Index of first position of current random vector
  
  for(i in 1:k){
    sumdim = sum(dim[1:i])
    eigen_Rii[1:dim[i],i] = pmax(eigen(R[start:sumdim,start:sumdim])$values,0)
    start = sumdim + 1
  }
  
  eigen_prod = pmax(eigen(Re(expm::sqrtm(R0)) %*% R %*% Re(expm::sqrtm(R0)))$values,0)
  num = s - sum(sqrt(eigen_prod)) # Numerator of D2
  denom = s - sum(sqrt(rowSums(eigen_Rii^2))) # Denominator of D2
  
  return(num/denom)
  
}
