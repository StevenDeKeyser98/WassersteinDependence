  # Returns the independence correlation matrix of the random vectors of dimensions dim corresponding to the correlation matrix R

create_R0 = function(R,dim){
  
  R0 = matrix(0,nrow(R),ncol(R))
  start = 1 # Index of first position of current random vector
  
  for(i in 1:length(dim)){
    sumdim = sum(dim[1:i]) # Index of last position of current random vector
    R0[start:sumdim,start:sumdim] = R[start:sumdim,start:sumdim]
    start = sumdim + 1
  }
  
  return(R0)
    
}
