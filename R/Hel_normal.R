# Function to compute Gaussian copula Hellinger distance given correlation matrix R
# with dimensions of the random vectors specified in dim

Hel_normal = function(R,dim){
  
  s = nrow(R)
  B = solve(create_R0(R,dim))
  det = det(diag(s) + B %*% R)
  
  return(1 - (sqrt((2^s)/det) * exp((-1/2)*MI_normal(R,dim))))
  
}
