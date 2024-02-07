# Function to compute Gaussian copula mutual information given correlation matrix R
# with dimensions of the random vectors specified in dim

MI_normal = function(R,dim){
  
  R0 = create_R0(R,dim)
  
  return((-1/2) * log(det(R)/det(R0)))
  
}
