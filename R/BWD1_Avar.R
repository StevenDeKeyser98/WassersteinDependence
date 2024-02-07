# Asymptotic standard deviation of first BW coefficient D1
# R is the correlation matrix 
# Dim should be in ascending order

BWD1_Avar = function(R,dim){
  
  s = nrow(R) # Total dimension
  k = length(dim) # Amount of random vectors
  eigen_Rii = matrix(0,s,k) # Eigenvalues of diagonal blocks
  start = 1 # Index of first position of current random vector
  L_matrices = list() # Lambda matrices
  U_matrices = list() # U matrices
  D1_matrices = list() # Delta matrices
  
  # Compute eigenvalues (Lambda matrices) and eigenvectors (U matrices)
  
  for(i in 1:k){ # Iterate over all random vectors 
    sumdim = sum(dim[1:i]) # Index of last position of current random vector
    eigen = eigen(R[start:sumdim,start:sumdim])
    L_matrices[[paste("L", i, i, sep = "")]] = diag(pmax(eigen$values,0))
    U_matrices[[paste("U", i, i, sep = "")]] = eigen$vectors
    eigen_Rii[1:dim[i],i] = pmax(eigen$values,0)
    start = sumdim + 1
  }
  
  # Compute Delta matrices
  
  sum = 0
  
  for(j in 1:k){sum = sum + L_matrices[[paste("L", j, j, sep = "")]][1:dim[1],1:dim[1]]}
  D1_matrices[["D1_1"]] = Re(expm::sqrtm(solve(sum))) # Delta1
  
  for(i in 2:k){ # Deltai for i = 2,...,k
    
    if(dim[i] == dim[i-1]){
      D1_matrices[[paste("D1_", i, sep = "")]] = D1_matrices[[paste("D1_", i-1, sep = "")]]
    } else{
      sum = 0
      
      for(j in i:k){sum = sum + L_matrices[[paste("L", j, j, sep = "")]][(dim[i-1]+1):dim[i],(dim[i-1]+1):dim[i]]}
      D1_matrices[[paste("D1_", i, sep = "")]] = adiag(D1_matrices[[paste("D1_", i-1, sep = "")]],
                                                       Re(expm::sqrtm(solve(sum))))
    }
  }
  
  R0 = create_R0(R,dim) # R0
  sqrtR0 = Re(expm::sqrtm(R0)) # Square root of R0
  invsqrtR0 = solve(sqrtR0) # Inverse square root of R0
  
  C1 = sum(sqrt(eigen_Rii)) - sum(sqrt(rowSums(eigen_Rii))) # C1
  Ups1 = U_matrices[["U11"]] %*% D1_matrices[["D1_1"]] %*% t(U_matrices[["U11"]]) # First block of Upsilon1
  
  for(i in 2:k){ # Construct other blocks of Upsilon1
    Ups1 = adiag(Ups1,U_matrices[[paste("U", i, i, sep = "")]] %*% D1_matrices[[paste("D1_", i, sep = "")]] %*% t(U_matrices[[paste("U", i, i, sep = "")]]))
  }
  
  D1 = BWD1(R,dim)
  M1 = (1/(2*C1)) * (-Re(expm::sqrtm(solve(R))) + (1-D1) * invsqrtR0 + D1 * Ups1)
  zeta1H = R%*%(M1 - diag(diag(M1 %*% R)))
  zeta1 = sqrt(2 * sum(diag(zeta1H %*% zeta1H)))
  
  return(zeta1)
  
}
