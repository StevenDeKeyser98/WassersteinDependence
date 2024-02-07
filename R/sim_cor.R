# Generate a correlation matrix with zero blocks

sim_cor = function(m,p,k){
  
  nmbrs = integer(p^2) # Must have length p^2
  
  for(i in 0:(p-1)){
    nmbrs[i*p + (i+1)] = 1
  }
  
  Q1 = kronecker(Matrix(nmbrs,p,p),diag(m))
  
  # Put upper diagonal elements equal to zero
  
  trilled = tril(Q1)
  
  # Add k colmuns of zeroes
  
  Q2 = cbind(trilled,Matrix(0,m*p,k))
  
  # Add k rows of length m*p + k 
  
  Q3 = rbind(Q2,cbind(Matrix(rnorm(k*m*p),k,m*p),Diagonal(k)))
  covmat = round(as.matrix(tcrossprod(Q3)),2)
  covmat = covmat + 0.05 * diag(rep(1, m*p + k))
  D = sqrt(solve(diag(diag(covmat))))
  cormat = D %*% covmat %*% D 
  
  return(cormat)
  
}
