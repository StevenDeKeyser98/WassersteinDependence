# Perform group lasso estimation 
# Sigma is an initial guess for the covariance matrix 
# S is ML estimator of the covariance matrix
# Model selection is done by choosing penalty parameter in omega such that BIC is maximal
# Dim contains the dimensions of the random vectors
# Step.size is the step size used in generalized gradient descent, affects the speed of the algorithm

GroupLasso = function(Sigma, S, n, omegas, dim, step.size){
  
  q = nrow(S) # Total dimensions
  bic = integer(length(omegas)) # BIC values
  ests = list() # Estimates
  
  for(o in 1:length(omegas)){
    omega = omegas[o]
    gspcov = Gspcov(Sigma, S, omega, dim, step.size) # Perform group lasso estimation
    names(gspcov)[[2]] = "sigma"
    sigma = gspcov$sigma
    penal = df(sigma,S,dim) * log(n) # Penalty
    bic[o] = -n*(log(det(sigma)) + sum(diag(solve(sigma) %*% S))) - penal # BIC
    ests[[o]] = gspcov
  }
  best = which.max(bic)
  
  return(list(ests[[best]],omegas[best]))
  
}

df = function(Sigma,S,dim){
  
  # Degrees of freedom for group lasso 
  # Sigma is candidate covariance matrix estimate
  # S is ML estimate of the covariance matrix
  # Dim contains the dimensions of the random vectors
  # See page 28 of paper 
  
  q = nrow(S) # Total dimension
  df = 0
  start = 1 # Index of first position of current random vector
  
  # Upper off-diagonal blocks 
  for(i in 1:(length(dim)-1)){
    sumdim = sum(dim[1:i]) # Index of last position of current random vector
    
    for(j in i:(length(dim)-1)){
      sumdim2 = sum(dim[1:j]) # Index of first position of next random vector - 1
      blockSigma = Sigma[start:sumdim,(sumdim2 + 1):(sumdim2 + dim[j+1])] # Sigma block
      blockS = S[start:sumdim,(sumdim2 + 1):(sumdim2 + dim[j+1])] # S block
      normblockSigma = sqrt(sum(blockSigma^2)) # Frobenius norm of Sigma block
      normblockS = sqrt(sum(blockS^2)) # Frobenius norm of S block
      df = df + as.numeric(normblockSigma > 0) * (1 + ((normblockSigma/normblockS) * ((dim[i] * dim[j+1]) - 1)))
    }
    start = sumdim + 1
  }
  start = 1
  
  # Diagonal blocks
  for(i in 1:length(dim)){
    sumdim = sum(dim[1:i])
    
    if(dim[i] == 1){
      df = df
    } else{
      delta = matrix(1,dim[i],dim[i]) - diag(dim[i]) # Deltai matrix
      blockSigma = Sigma[start:sumdim,start:sumdim]
      blockS = S[start:sumdim,start:sumdim]
      normblockSigma = sqrt(sum((delta * blockSigma)^2))
      normblockS = sqrt(sum((delta * blockS)^2))
      df = df + as.numeric(normblockSigma > 0) * (1 + ((normblockSigma/normblockS) * (((dim[i] * (dim[i] - 1))/2) - 1)))
    }
    start = sumdim + 1
  }
  df = df + q
  
  return(df)
}

Gspcov = function(Sigma, S, omega, dim, step.size, nesterov = TRUE,
                  n.outer.steps = 1e4, n.inner.steps = 1e4,
                  tol.outer = 1e-4, thr.inner = 1e-2,
                  backtracking = 0.2, trace = 0){
  
  # Computes the group lasso estimator where Frobenius penalty is imposed 
  # on all off-diagonal blocks, and all diagonal blocks, but no penalization of
  # diagonal elements
  # Sigma is an initial guess for the covariance matrix 
  # S is ML estimator of the covariance matrix
  # Dim contains the dimensions of the random vectors
  # See R package spcov for explanations of other arguments and code
  
  if(all(omega == 0)){
    cat("Skipping MM.  Solution is S!", fill = T)
    
    return(list(n.iter=0, Sigma=S, obj=ComputeObjective(S, S, omega, dim = dim)))
  }
  
  stopifnot(omega >= 0)
  
  if(trace > 0){
    cat("---", fill = T)
    cat(ifelse(nesterov, "using Nesterov, ", ""))
    cat(ifelse(backtracking, "backtracking line search", ""), fill = T)
    cat("---", fill = T)
  }
  
  mean.abs.S = mean(abs(S))
  
  if(min(eigen(Sigma, symmetric = T, only.values = T)$val) < 1e-5)
    
    warning("Starting value is nearly singular.")
  
  del = ComputeDelta(S, omega, trace = trace - 1, dim = dim) 
  objective = ComputeObjective(Sigma, S, omega, dim = dim)
  
  if(trace > 0)
    
    cat("objective: ", objective, fill = T)
  
  n.iter = NULL 
  
  for(i in seq(n.outer.steps)){
    Sigma0 = Sigma
    
    if(trace > 0)
      cat("step size given to GGDescent/Nesterov:", step.size, fill = T)
    gg = GGDescent(Sigma = Sigma, Sigma0 = Sigma0, S = S, omega = omega,
                   del = del, nsteps = n.inner.steps,
                   step.size = step.size,
                   nesterov = nesterov,
                   tol = thr.inner * mean.abs.S,
                   trace = trace - 1,
                   backtracking = backtracking, dim = dim)
    
    Sigma = gg$Sigma
    objective = c(objective, ComputeObjective(Sigma, S, omega, dim = dim))
    
    if(trace > 0){
      cat("objective: ", objective[length(objective)],
          " (", gg$niter, "iterations, max step size:",
          max(gg$step.sizes), ")", fill = T)
    }
    
    if(backtracking){
      
      if(max(gg$step.sizes) < step.size * backtracking ^ 2){
        step.size = step.size * backtracking
        
        if(trace > 0)
          cat("Reducing step size to", step.size, fill = T)
      }
    }
    n.iter = c(n.iter, gg$niter)
    
    if(objective[i + 1] > objective[i] - tol.outer){
      cat("MM converged in", i, "steps!", fill = T)
      break
    }
  }
  
  list(n.iter = n.iter, Sigma = gg$Sigma, obj = objective)
}

GGDescent = function(Sigma, Sigma0, S, omega, del, nsteps,
                     step.size, nesterov = FALSE, backtracking = FALSE,
                     tol = 1e-3, trace = 0, dim){
  
  # Generalized gradient descent steps
  # See spcov package for explanation of most of the code
  # For group lasso, the ComputeObjective function is different
  
  if(backtracking){
    beta = backtracking
    
    if (beta <= 0 | beta >= 1)
      stop("Backtracking parameter beta must be in (0,1).")
  }
  tt = step.size
  converged = FALSE
  exit = FALSE
  obj.starting = ComputeObjective(Sigma, S, omega, dim = dim)
  Sigma.starting = Sigma
  Omega = Sigma
  Sigma.last = Sigma
  ttts = NULL
  ttt = tt 
  
  for(i in seq(nsteps)){
    inv.Sigma0 = solve(Sigma0)
    log.det.Sigma0 = LogDet(Sigma0)
    grad.g = ComputeGradientOfg(Omega, S, Sigma0, inv.Sigma0 = inv.Sigma0)
    grad.g = (grad.g + t(grad.g)) / 2 
    g.omega = g(Omega, S, Sigma0, inv.Sigma0 = inv.Sigma0, log.det.Sigma0 = log.det.Sigma0)
    
    while(backtracking){
      soft.thresh = ProxADMM(Omega - ttt * grad.g, del, omega = omega * ttt, dim = dim, rho = .1)$X
      gen.grad.g = (Omega - soft.thresh) / ttt
      left = g(soft.thresh, S, Sigma0,inv.Sigma0 = inv.Sigma0, log.det.Sigma0 = log.det.Sigma0)
      right = g.omega - ttt * sum(grad.g * gen.grad.g) + ttt * sum(gen.grad.g ^ 2) / 2
      
      if(is.na(left) || is.na(right)){
        print("left or right is NA.")
        # browser()
      }
      
      if(left <= right){
        Sigma = soft.thresh
        ttts = c(ttts, ttt)
        
        if(mean(abs(Sigma - Sigma.last)) < tol){
          converged = TRUE
          break 
        }
        
        if(nesterov)
          Omega = Sigma + (i - 1) / (i + 2) * (Sigma - Sigma.last)
        
        else
          Omega = Sigma
        Sigma.last = Sigma
        
        if(trace > 0)
          cat("--true objective:", ComputeObjective(Sigma, S, omega, dim = dim), fill = T)
        
        if(trace > 0)
          cat(i, ttt, " ")
        break
      }
      ttt = beta * ttt
      
      if(ttt < 1e-15){
        cat("Step size too small: no step taken", fill = T)
        exit = TRUE
        break
      }
    }
    
    if(!backtracking){
      Sigma = ProxADMM(Sigma - ttt * grad.g, del, omega = omega * ttt, dim = dim, rho=.1)$X
      
      if(mean(abs(Sigma - Sigma.last)) < tol)
        converged = TRUE
      
      if(nesterov)
        Omega = Sigma + (i - 1)/(i + 2) * (Sigma - Sigma.last)
      
      else
        Omega = Sigma
      Sigma.last = Sigma
    }
    
    if(converged){
      
      if(trace > 0){
        cat("--GG converged in", i, "steps!")
        
        if(backtracking)
          cat(" (last step size:", ttt, ")", fill = T)
        
        else
          cat(fill = T)
      }
      break
    }
    
    if(exit){
      break
    }
  }
  obj.end = ComputeObjective(Sigma, S, omega, dim = dim)
  
  if(obj.starting < obj.end){
    
    if(nesterov){
      cat("Objective rose with Nesterov.  Using generalized gradient instead.", fill = T)
      
      return(GGDescent(Sigma = Sigma.starting, Sigma0 = Sigma0, S = S, omega = omega,
                       del = del, nsteps = nsteps, step.size = step.size,
                       nesterov = FALSE,
                       backtracking = backtracking,
                       tol = tol, trace = trace, dim = dim))
    }
    # browser()
    cat("--Returning initial Sigma since GGDescent/Nesterov did not decrease objective", fill=T)
    Sigma = Sigma.starting
  }
  list(Sigma = Sigma, niter = i, step.sizes = ttts)
}


ComputeObjective = function(Sigma, S, omega, dim){
  
  # Group lasso objective function
  
  -2 * ComputeLikelihood(Sigma, S) + ComputePenalty(Sigma, omega, dim = dim)
}

ComputeLikelihood = function(Sigma, S){
  
  # Gaussian likelihood
  
  p = nrow(Sigma)
  ind.diag = 1 + 0:(p - 1) * (p + 1)
  -(1 / 2) * (LogDet(Sigma) + sum(solve(Sigma, S)[ind.diag]))
}

ComputePenalty = function(Sigma, omega, dim){
  
  # Group lasso penalty
  
  penal = 0
  start = 1 # Index of first position of current random vector
  
  # Upper off-diagonal blocks 
  for(i in 1:(length(dim)-1)){
    sumdim = sum(dim[1:i]) # Index of last position of current random vector
    
    for(j in i:(length(dim)-1)){
      sumdim2 = sum(dim[1:j]) # Index of first position of next random vector - 1
      block = Sigma[start:sumdim,(sumdim2 + 1):(sumdim2 + dim[j+1])] # Sigma block
      penal = penal + sqrt(dim[i] * dim[j+1] * sum(block^2))
    }
    start = sumdim + 1
  }
  penal = 2 * omega * penal
  start = 1
  penal2 = 0
  
  # Diagonal blocks
  for(i in 1:length(dim)){
    sumdim = sum(dim[1:i])
    delta = matrix(1,dim[i],dim[i]) - diag(dim[i])
    block = Sigma[start:sumdim,start:sumdim]
    penal2 = penal2 + sqrt(dim[i] * (dim[i] - 1) * sum((delta * block)^2))
    start = sumdim + 1
  }
  
  return(penal + omega * penal2)
}

LogDet = function(M) {
  
  # Log of determinant of matrix M, see spcov package
  
  log.det = determinant(M, logarithm=TRUE)
  
  if(log.det$sign == -1)
    return(NA)
  as.numeric(log.det$mod)
}

h = function(x, a) log(x) + a / x

dhdx = function(x, a) (1 - a / x) / x

FindMinSolution = function(a, c, tol = 1e-10, max.iter = 1e3, trace = 1){
  
  # See spcov package
  
  if (h(a, a) > c)
    stop("No solution!")
  
  if (h(a, a) == c)
    return(a)
  
  
  x.min = 0.1 * a
  
  for(i in seq(max.iter)){
    
    if(trace > 1)
      cat("h(", x.min, ", a) = ", h(x.min, a), fill=T)
    x.min = x.min - (h(x.min, a) - c) / dhdx(x.min, a)
    
    if(x.min <= 0)
      x.min = a * tol
    
    else if(x.min >= a)
      x.min = a * (1 - tol)
    
    if(abs(h(x.min, a) - c) < tol & h(x.min, a) >= c){
      
      if(trace > 0)
        cat("Eval-bound converged in ", i, "steps to ", x.min, fill=T)
      break
    }
  }
  x.min
}

ComputeDelta = function(S, omega, Sigma.tilde = NULL, trace = 1, dim){
  
  # See spcov package
  
  if(is.null(Sigma.tilde))
    Sigma.tilde = diag(diag(S))
  
  p = nrow(S)
  f.tilde = ComputeObjective(Sigma.tilde, S, omega, dim = dim)
  minev = min(eigen(S)$val)
  c = f.tilde - (p - 1) * (log(minev) + 1)
  FindMinSolution(a = minev, c = c, trace = trace)
}

ComputeGradientOfg = function(Sigma, S, Sigma0, inv.Sigma0=NULL){
  
  # See spcov package
  
  inv.Sigma = solve(Sigma)
  
  if(is.null(inv.Sigma0))
    solve(Sigma0) - inv.Sigma %*% S %*% inv.Sigma
  
  else
    inv.Sigma0 - inv.Sigma %*% S %*% inv.Sigma
}

g = function(Sigma, S, Sigma0, inv.Sigma0 = NULL, log.det.Sigma0 = NULL) {
  
  # See spcov package
  
  p = nrow(Sigma)
  ind.diag = 1 + 0:(p - 1) * (p + 1)
  
  if (is.null(log.det.Sigma0))
    log.det.Sigma0 = LogDet(Sigma0)
  
  if (is.null(inv.Sigma0))
    log.det.Sigma0 + sum((solve(Sigma0, Sigma) + solve(Sigma, S))[ind.diag]) - p
  
  else
    log.det.Sigma0 + sum((inv.Sigma0 %*% Sigma + solve(Sigma, S))[ind.diag]) - p
}

ProxADMM = function(A, del, omega, dim, rho = .1, tol = 1e-6, maxiters = 100, verb = FALSE){
  
  # See spcov package
  
  soft = SoftThreshold(A, omega, dim)
  minev = min(eigen(soft, symmetric = T, only.values = T)$val)
  
  if(minev >= del){
    return(list(X = soft, Z = soft, obj = ComputeProxObjective(soft, A, omega, dim = dim)))
  }
  
  p = nrow(A)
  obj = NULL
  
  Z = soft
  Y = matrix(0, p, p)
  
  for(i in seq(maxiters)){    
    B = (A + rho * Z - Y) / (1 + rho)
    
    if(min(eigen(B, symmetric = T, only.values = T)$val) < del){
      eig = eigen(B, symmetric = T)
      X = eig$vec %*% diag(pmax(eig$val, del)) %*% t(eig$vec)
    }
    
    else{
      X = B
    }
    obj = c(obj, ComputeProxObjective(X, A, omega, dim = dim))
    
    if(verb)
      cat(" ", obj[i], fill = T)
    
    if(i > 1)
      
      if(obj[i] > obj[i - 1] - tol){
        
        if(verb)
          cat(" ADMM converged after ", i, " steps.", fill = T)
        break
      }
    
    Z = SoftThreshold(X + Y / rho, omega / rho, dim = dim)
    Y = Y + rho * (X - Z)
  }
  
  list(X = X, Z = Z, obj = obj)
}

SoftThreshold = function(x, omega, dim){
  
  # Elementwise soft thresholding operation in case of group lasso
  # X is ML estimator of covariance matrix
  # Omega is the penalty parameter
  # Dim contains the dimensions of the random vectors
  
  q = nrow(x) # Total dimension
  Sol = matrix(0,q,q) # Solution 
  start = 1 # Index of first position of current random vector
  
  # Upper off-diagonal blocks 
  for(i in 1:(length(dim)-1)){
    sumdim = sum(dim[1:i]) # Index of last position of current random vector
    
    for(j in i:(length(dim)-1)){
      sumdim2 = sum(dim[1:j]) # Index of first position of next random vector - 1
      block = x[start:sumdim,(sumdim2 + 1):(sumdim2 + dim[j+1])]
      S = sqrt(sum(block^2)) ; gamma = sqrt(dim[i] * dim[j+1]) # Si,m and gammai,m
      
      if(S <= omega * gamma){
        Sol[start:sumdim,(sumdim2 + 1):(sumdim2 + dim[j+1])] = matrix(0,dim[i],dim[j+1])
      } else{
        Sol[start:sumdim,(sumdim2 + 1):(sumdim2 + dim[j+1])] = block * (1 - ((omega * gamma) / S))
      }
      Sol[(sumdim2 + 1):(sumdim2 + dim[j+1]),start:sumdim] = t(Sol[start:sumdim,(sumdim2 + 1):(sumdim2 + dim[j+1])])
    }
    start = sumdim + 1
  }
  start = 1
  
  # Diagonal blocks
  for(i in 1:length(dim)){
    sumdim = sum(dim[1:i])
    delta = matrix(1,dim[i],dim[i]) - diag(dim[i])
    block = x[start:sumdim,start:sumdim]
    S = sqrt(sum((delta * block)^2)) ; gamma = sqrt(dim[i] * (dim[i] - 1))
    
    if(S <= omega * gamma){
      Sol[start:sumdim,start:sumdim] = matrix(0,dim[i],dim[i])
    } else{
      Sol[start:sumdim,start:sumdim] = block * (1 - ((omega * gamma) / S))
    }
    
    if(length(start:sumdim) == 1){
      Sol[start:sumdim,start:sumdim] = block
    } else{
      diag(Sol[start:sumdim,start:sumdim]) = diag(block)
    }
    start = sumdim + 1
  }

  return(Sol)

}

ComputeProxObjective = function(X, A, omega, dim){
  
  # Objective (29) of paper
  
  sum((X-A)^2) / 2 + ComputePenalty(X, omega, dim)

}
