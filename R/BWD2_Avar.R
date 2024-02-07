# Asymptotic standard deviation of second BW coefficient D2
# R is the correlation matrix 
# Dim should be in ascending order

BWD2_Avar = function(R,dim){
  
  s = nrow(R) # Total dimension
  k = length(dim) # Amount of random vectors
  eigen_Rii = matrix(0,s,k) # Eigenvalues of diagonal blocks
  start = 1 # Index of first position of current random vector
  L_matrices = list() # Lambda matrices
  U_matrices = list() # U matrices
  D2_matrices = list() # Delta tilde matrices
  
  # Compute eigenvalues (Lambda matrices) and eigenvectors (U matrices)
  
  for(i in 1:k){ # Iterate over all random vectors
    sumdim = sum(dim[1:i]) # Index of last position of current random vector 
    eigen = eigen(R[start:sumdim,start:sumdim])
    L_matrices[[paste("L", i, i, sep = "")]] = diag(pmax(eigen$values,0))
    U_matrices[[paste("U", i, i, sep = "")]] = eigen$vectors
    eigen_Rii[1:dim[i],i] = pmax(eigen$values,0)
    start = sumdim + 1
  }
  
  # Compute Delta tilde matrices
  
  sum = 0
  
  for(j in 1:k){sum = sum + L_matrices[[paste("L", j, j, sep = "")]][1:dim[1],1:dim[1]]^2}
  D2_matrices[["D2_1"]] = Re(expm::sqrtm(solve(sum))) %*% L_matrices[["L11"]][1:dim[1],1:dim[1]] # Delta tilde1
  
  for(i in 2:k){ # Delta tildei for i = 2,...,k
    D = solve(L_matrices[[paste("L", i-1, i-1, sep = "")]][1:dim[1],1:dim[1]]) %*% L_matrices[[paste("L", i, i, sep = "")]][1:dim[1],1:dim[1]]
    
    if(i > 2){
      
      for(j in 2:(i-1)){
        
        if(dim[j] != dim[j-1]){
          D = adiag(D,solve(L_matrices[[paste("L", i-1, i-1, sep = "")]][(dim[j-1]+1):dim[j],(dim[j-1]+1):dim[j]]) %*% 
                    L_matrices[[paste("L", i, i, sep = "")]][(dim[j-1]+1):dim[j],(dim[j-1]+1):dim[j]])
        }
      }
    }
    
    if(dim[i] == dim[i-1]){
      D2_matrices[[paste("D2_", i, sep = "")]] = D2_matrices[[paste("D2_", i-1, sep = "")]] %*% D
    } else{
      sum = 0
      
      for(j in i:k){sum = sum + L_matrices[[paste("L", j, j, sep = "")]][(dim[i-1]+1):dim[i],(dim[i-1]+1):dim[i]]^2}
      D2_matrices[[paste("D2_", i, sep = "")]] = adiag(D2_matrices[[paste("D2_", i-1, sep = "")]] %*% D,
                                                       Re(expm::sqrtm(solve(sum))) %*% 
                                                       L_matrices[[paste("L", i, i, sep = "")]][(dim[i-1]+1):dim[i],(dim[i-1]+1):dim[i]])
    }
  }
  
  R0 = create_R0(R,dim) # R0
  sqrtR0 = Re(expm::sqrtm(R0)) # Square root of R0
  invsqrtR0 = solve(sqrtR0) # Inverse square root of R0
  J = invsqrtR0 %*% Re(expm::sqrtm(sqrtR0 %*% R %*% sqrtR0)) %*% invsqrtR0 # J
  J0 = create_R0(J,dim) # J0

  C2 = s - sum(sqrt(rowSums(eigen_Rii^2))) # C2
  Ups2 = U_matrices[["U11"]] %*% D2_matrices[["D2_1"]] %*% t(U_matrices[["U11"]]) # First block of Upsilon2
  
  for(i in 2:k){ # Construct other blocks of Upsilon2
    Ups2 = adiag(Ups2,U_matrices[[paste("U", i, i, sep = "")]] %*% D2_matrices[[paste("D2_", i, sep = "")]] %*% t(U_matrices[[paste("U", i, i, sep = "")]]))
  }
  
  D2 = BWD2(R,dim)
  M2 = (1/C2) * ((-1/2) * (J0 + solve(J)) + (1-D2) * diag(rep(1,s)) + D2 * Ups2)
  zeta2H = R%*%(M2 - diag(diag(M2 %*% R)))
  zeta2 = sqrt(2 * sum(diag(zeta2H %*% zeta2H)))
  
  return(zeta2)
  
}
