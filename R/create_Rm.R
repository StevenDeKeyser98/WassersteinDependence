# Returns the correlation matrix with maximal BW dependence of the random vectors
# of dimensions dim corresponding to the independence correlation matrix R0

create_Rm = function(R0,dim){
  
  Rm = matrix(0,nrow(R0),ncol(R0))
  start = 1 # Index of first position of current random vector
  
  for(i in 1:length(dim)){
    sumdim = sum(dim[1:i]) # Index of last position of current random vector
    Rm[start:sumdim,start:sumdim] = R0[start:sumdim,start:sumdim] # Diagonal blocks
    
    start2 = sumdim + 1 # Index of first position of next random vector 
      
    if(i < length(dim)){
        
      for(j in (i+1):length(dim)){ # Only look at upper diagonal blocks
        sumdim2 = sum(dim[1:j]) # Index of last position of next random vector
        block = R0[start:sumdim,start:sumdim] # Diagonal block of first random vector
        block2 = R0[start2:sumdim2,start2:sumdim2] # Diagonal block of second random vector
        ident = diag(length(start:sumdim) + length(start2:sumdim2))
        ident = ident[1:length(start:sumdim), 1:length(start2:sumdim2)] # Piij
        spectral = eigen(block)
        spectral2 = eigen(block2)
        Psiblock = spectral$vectors %*% diag(sqrt(spectral$values)) %*%
          ident %*% diag(sqrt(spectral2$values)) %*% t(spectral2$vectors)
          
        Rm[start:sumdim,start2:sumdim2] = Psiblock # Psiij
        Rm[start2:sumdim2,start:sumdim] = t(Rm[start:sumdim,start2:sumdim2])
        start2 = sumdim2 + 1
      }
    }
    start = sumdim + 1
  }
  
  return(Rm)
      
}
