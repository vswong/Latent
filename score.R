#' Gradient of the log-likelihood
#' 
#' Evaluates the score functions (gradients of the log-likelihood) for a 
#' Poisson distribution with log-linear model for the mean.
#'
#' @param data matrix of observations, with rows equal to the number of 
#'     observations and columns equal to the number of categories
#' @param params vector of model parameters, consisting of 
#'     \eqn{\boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\gamma}}.
#' @param event vector of labels with one entry per row of data, where 
#'     each entry identifies the event to which the corresponding row belongs
#' @param specific vector of logical values, one entry per column of data, 
#'     with each entry indicating whether the corresponding column of the 
#'     data is a human-specific bacterium
#' 
#' @return The score functions, evaluated at the current parameter values.
#' 
score <- function(data, params, event, specific=NULL) {
  # basic constants:
  n <- nrow(data)
  p <- ncol(data)
  d <- length(unique(event))
  
  # check for predictable error conditions
  if (is.null(specific))
    specific <- rep(FALSE, p)
  if (length(specific) != p)
    stop("Argument 'specific' must have length equal to the number of 
         pathogens.")
  if (!all(typeof(specific) == "logical"))
    stop("Each element of 'specific' must be logical (boolean).")
  
  # Extract parameters from the parameter vector.
  # sigma is length 2
  # lambda is the lagrange multiplier
  alpha <- params[1:(p * d)]
  beta <- params[(p * d + 1):(p * (2 + d))]
  gamma <- params[(p * (2 + d) + 1):(p * (2 + d) + 2 * n)]
  lambda <- tail(params, d)

  # read the separate elements of beta and gamma
  # beta is 2 per column
  # gamma is 2 per row
  beta1 <- beta[1:p]
  beta2 <- tail(beta, p)
  gamma1 <- gamma[1:n]
  gamma2 <- tail(gamma, n)

  
  # event-level alphas
  # sort out lambdas here too
  alpha.local <- matrix(0, n, p)
  lambda.local <- vector()
  gamma1.partial <- gamma1
  gamma2.partial <- gamma2
  for (k in 1:d) {
    indx <- which(event==unique(event)[k])
    alpha.local[indx,] <- matrix(alpha[(p * (k - 1) + 1):(p * k)], length(indx), p, byrow=TRUE)
    
    #Derivative w.r.t. lambda(lagrange multiplier) is simply the dot product of the gammas in the event
    lambda.local <- c(lambda.local, sum(gamma1[indx]*gamma2[indx])^2)
    
    gamma1.partial[indx] = -lambda[k]*gamma2[indx]
    gamma2.partial[indx] = -lambda[k]*gamma1[indx]
    if(sum(gamma1[indx]*gamma2[indx]) < 0){
      gamma1.partial[indx] = -gamma1.partial[indx]
      gamma2.partial[indx] = -gamma2.partial[indx]
    }
  }
  
  # compute this once to save time:
  eta <- exp(gamma1 %*% t(beta1) + gamma2 %*% t(beta2) + alpha.local)
  
  # gradient in the direction of event-specific alphas:
  grad <- vector()
  #Padding for alphas/betas
  grad <- c(grad, rep(0, d*p))
  grad <- c(grad, rep(0, p))
  grad <- c(grad, rep(0, p))
  
  # gradient in the direction of gamma:
  # cat(paste("Magnitude for log-likelihood:", sqrt(sum(sweep(data - eta, 2, beta1, '*'), na.rm=TRUE)^2), "Magnitude for lagrange:", sqrt(sum(gamma1.partial^2)), "\n"))
  grad <- c(grad, rowSums(sweep(data - eta, 2, beta1, '*'), na.rm=TRUE)+gamma1.partial)
  grad <- c(grad, rowSums(sweep(data - eta, 2, beta2, '*'), na.rm=TRUE)+gamma2.partial)
  
  
  # lambdas as well
  grad <- c(grad, c(0,0,0,0))
  
  grad
}
