# Computes the first BW dependence coefficient given the correlation matrix R
# of the random vectors of dimensions dim

BWD1 = function(R,dim){
  
  s = nrow(R) # Total dimension
  k = length(dim) # Amount of random vectors 
  eigen_R = pmax(eigen(R)$values,0) # Eigenvalues of R
  eigen_Rii = matrix(0,s,k) # # Eigenvalues of diagonal blocks
  start = 1 # Index of first position of current random vector
  
  for(i in 1:k){
    sumdim = sum(dim[1:i])
    eigen_Rii[1:dim[i],i] = pmax(eigen(R[start:sumdim,start:sumdim])$values,0) 
    start = sumdim + 1
  }
  
  num = sum(sqrt(eigen_Rii)) - sum(sqrt(eigen_R)) # Numerator of D1
  denom = sum(sqrt(eigen_Rii)) - sum(sqrt(rowSums(eigen_Rii))) # Denominator of D1
  
  return(num/denom)
  
}
