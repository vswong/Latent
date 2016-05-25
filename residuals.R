residuals <- function(result){
  n = length(result$data$mei)
  p = ncol(result$alpha)
  d = unique(result$event)
  
  alpha.local = matrix(0, n, p)
  for (k in d) {
    indx = which(result$event==unique(result$event)[k])
    alpha.local[indx,] = matrix(as.vector(t(result$alpha))[p*(k-1) + 1:p], length(indx), p, byrow=TRUE)
    
  }
  mu = round(exp(alpha.local+result$gamma%*%t(result$beta)))
  mu[mu==0] = 1
  pearsons = (round(result$data) - mu)/sqrt(mu)
  
  deviance = 2*(result$data*(log(result$data+1)-log(mu))-(result$data-mu))
  
  
  
  
  
  
  
  
  n = length(result2$data$mei)
  p = ncol(result2$alpha)
  d = unique(result2$event)
  
  alpha.local = matrix(0, n, p)
  for (k in d) {
    indx = which(result2$event==unique(result2$event)[k])
    alpha.local[indx,] = matrix(as.vector(t(result2$alpha))[p*(k-1) + 1:p], length(indx), p, byrow=TRUE)
    
  }
  mu = round(exp(alpha.local+result2$gamma1%*%t(result2$beta1)+result2$gamma2%*%t(result2$beta2)))
  mu[mu==0] = 1
  pearsons2 = (round(result2$data) - mu)/sqrt(mu)
  
  deviance2 = 2*(result2$data*(log(result2$data+1)-log(mu))-(result2$data-mu))
}