# Sample normal scores rank correlation matrix 
# Sample is a matrix with the observations in the rows
# Lambda is the Ridge penalty parameter (1 for no penalization)

est_R = function(sample,lambda){
  
  n = nrow(sample) # Sample size
  s = ncol(sample) # Total dimension 
  
  scores = matrix(0,n,s) # Normal scores 
  
  for(j in 1:s){
    scores[,j] = qnorm((n/(n+1)) * ecdf(sample[,j])(sample[,j]))
  }
  
  R_est = cor(scores)
  
  return(lambda * R_est + (1-lambda) * diag(s))
  
}
