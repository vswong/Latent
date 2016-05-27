irls <- function(data, params, event, specific) {
  #Basic constants:
  n <- nrow(data)
  p <- ncol(data)
  d <- length(unique(event))
  
  # Extract parameters from the parameter vector.
  # sigma is length 2
  # lambda is the lagrange multiplier
  alpha <- params[1:(p * d)]
  beta <- params[(p * d + 1):(p * (2 + d))]
  gamma <- params[(p * (2 + d) + 1):(p * (2 + d) + 2 * n)]
  lambda <- tail(params, d)
  
  # read the separate elements of beta and gamma
  # beta is 2 per column
  # gamma is 2 per row1111111
  beta1 <- beta[1:p]
  beta2 <- tail(beta, p)
  gamma1 <- gamma[1:n]
  gamma2 <- tail(gamma, n)
  
  #Make structure for alpha:
  alpha.local = matrix(0, n, p)
  
  #Make the design matrix here
  designMat.partial1 = matrix(0, n*p, p*d)
  designMat.partial2 = matrix(0, n*p, 2*p)
  for (k in 1:d) {
    indx <- which(event==unique(event)[k])
    alpha.local[indx,] = matrix(alpha[(p*(k-1)+1):(p*k)], length(indx), p, byrow=TRUE)
    for (j in 1:p){
      if(!specific[j])
      designMat.partial1[indx+(j-1)*n,1+(j-1)+(k-1)*p] = 1
    }
  }
  
#   designMat.partial2 = rbind(cbind(gamma1, matrix(0, n, 4), gamma2, matrix(0, n, 4)), 
#                              cbind(matrix(0,n,1), gamma1, matrix(0, n, 4), gamma2, matrix(0, n, 3)),
#                              cbind(matrix(0,n,2), gamma1, matrix(0, n, 4), gamma2, matrix(0, n, 2)),
#                              cbind(matrix(0,n,3), gamma1, matrix(0, n, 4), gamma2, matrix(0, n, 1)),
#                              cbind(matrix(0,n,4), gamma1, matrix(0, n, 4), gamma2))
  
  designMat.partial2 = rbind(cbind(gamma1, matrix(0, n, 4), gamma2, matrix(0, n, 4)), 
                             cbind(matrix(0,n,1), gamma1, matrix(0, n, 4), gamma2, matrix(0, n, 3)),
                             cbind(matrix(0,n,2), gamma1, matrix(0, n, 4), gamma2, matrix(0, n, 2)),
                             cbind(matrix(0,n,3), gamma1, matrix(0, n, 4), matrix(0, n, 1), matrix(0, n, 1)),
                             cbind(matrix(0,n,4), gamma1, matrix(0, n, 4), matrix(0, n, 1)))
  
  designMat = cbind(designMat.partial1, designMat.partial2)
  
  #Make eta into vector, should be 1x(5*38) = 1x190
  eta = as.vector(gamma1 %*% t(beta1) + gamma2 %*% t(beta2) + alpha.local)
  mu.local = exp(eta)
  wt = mu.local
  
  #Deduct alpha.local here as the design matrix will not include it
  z = eta + (as.vector(as.matrix(data)) - mu.local) / mu.local
  
  #Create response vector, consisting of vector z, 2 vectors of n length for the gaussian priors, and 1 entry for the lagrange multiplier
  response = as.vector(as.matrix(data))
  
  #Create initial "guess", consisting of alpha, beta1 and beta2
  start = c(alpha, beta1, beta2)
  
  pois = glm.fit(x=designMat, y=response, intercept=FALSE, family=poisson())
  
  alphabeta.new = pois$coefficients
  alphabeta.new[is.na(alphabeta.new)] = 0
  alphabeta.new
}