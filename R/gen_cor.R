# Generate random sparse correlation matrix 

gen_cor = function(q,sparsity){

  mask = rand(q)
  mask = mask * (mask < sparsity)
  mask[lower.tri(mask, diag = TRUE)] = 0
  mask = mask + t(mask) + eye(q)
  mask[mask > 0] = 1
  sigma = matrix(rnorm(q^2), q)
  sigma[lower.tri(sigma, diag = TRUE)] = 0
  sigma = sigma + t(sigma) + eye(q)
  sigma = sigma * mask
  sigma = sigma - (min(eig(sigma))-.1) * eye(q)
  D = sqrt(solve(diag(diag(sigma))))
  R = D %*% sigma %*% D 
  
  return(R)
  
}
