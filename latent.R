#' Estimate parameters of a latent variable model
#' 
#' Estimate the parameters of a latent variable model. Estimation is by the
#' method of maximum likelihood, with conjugate gradient descent used to 
#' maximize likelihood. An EM algorithm is used to impute the values for
#' censored observations.
#' 
#' @param data matrix of observed counts, number of rows equal to the number of
#'     observations and columns equal to the number of categories.
#' @param min.detect minimum limit of detection for the FIB assay - values
#'     below this threshold are censored
#' @param event vector of event assignment labels, one entry per row of data.
#'     Indicates the event to which each row of data belongs.
#' @param specific vector of TRUE/FALSE values, one entry per column of data.
#'     Each entry indicates whether the corresponding FIB species is
#'     human-specific.
#' @param verbose should the function provide verbose output about its
#'     progress? Defaults to TRUE.
#' 
#' @return A list containing the following elements:
#'     \item{data}{final version of the data, where censored values have been
#'         imputed by the EM algorithm (see \code{\link{E.step}})}
#'     \item{min.detect}{a vector where each element is the minimum level of
#'         detection for the corresponding column of data (copied from function
#'         input)}
#'     \item{event}{vector listing the event ID associated with the
#'         corresponding row of the data (copied from the input)}
#'     \item{specific}{vector indicating whether the corresponding column of
#'         the data is a human-specific indicator bacteria - in which case its
#'         intercept is set to zero (copied from input)}
#'     \item{alpha}{matrix of estimated intercept parameters for each indicator
#'         bacteria for each event - each row is for one event, and each column
#'         is for one species of bacteria}
#'     \item{beta}{vector of estimated slope parameters for the indicator
#'         bacteria, representing the expected marginal increase in log-
#'         concentration when the contamination index increases by one unit}
#'     \item{gamma}{vector of estimated contamination indices, one for each row
#'         of the data}
#' 
#' @examples # This script is an example of using the latent package to 
#' # estimate the parameters of a latent variable model for the
#' # FIB counts in the storm sewer data set (with event 3 removed).
#' 
#' # Import data, drop event 3, change mei[4] from "TNTC" to 0, and then
#' # convert all FIB to numerics:
#' data("dfOptAnalysisDataSSJan2015")
#' indx = which(dfOptAnalysisDataSSJan2015$Event != "03")
#' fib = dfOptAnalysisDataSSJan2015[indx, c("mei", "modmtec", "FC",
#'     "Bac.human", "Lachno.2")]
#' fib$mei[4] = 0
#' for (n in names(fib))
#'     fib[[n]] = as.numeric(fib[[n]])
#' 
#' # Set the censoring values for each of the FIB (these are guesstimates)
#' min.detect = c('mei'=1, 'modmtec'=1, 'FC'=1, 'Bac.human'=225,
#'     'Lachno.2'=225)
#' 
#' # The human-specific FIB are Bac.human and Lachno.2, which are the fourth
#' # and fifth columns
#' specific = c(FALSE, FALSE, FALSE, TRUE, TRUE)
#' 
#' # Get the event IDs
#' event = as.integer(dfOptAnalysisDataSSJan2015$Event[indx])
#' 
#' # Now estimate the model parameters:
#' latent(fib, min.detect, event, specific)
#' 
#' @export
latent <- function(data, min.detect, event, specific=NULL, verbose=TRUE) {
  #Initial parameters:
  n = nrow(data)
  p = ncol(data)
  d = length(unique(event))
  #alpha, beta, gamma, and sigma respectively
  #Mu is fixed at zero
  xx = c(rep(as.integer(!specific), length(unique(event))), rep(1, p), rep(1, n), 1)
  # xx = c(7.837, 2.667, 5.251, 0,0,8.002,3.77,7.25,0,0,7.609,3.208,6.314,0,0,1.316,3.612,2.574,5.519,5.227,1.905,1.888,1.934,1.749,2.001,2.048,2.015,1.781,0.972,1.055,1.751,1.726,1.657,1.453,1.639,1.009,1.114,0.954,1.013,0.782,0.597,1.2009,0.832,0.658,1.009,0.784,1.42,0.987,1.234,0.957,1.58,1.521,1.65,1.005,1.971,2.11,1.692,1.748,0.989)
  finished = FALSE
  
  tol = sqrt(.Machine$double.eps)
  tol = 1e-5
  check=Inf
  
  while (!finished) {
    #Loop for M-step:
    converged = FALSE
    f.mold = log.lik(data, xx, event)
    f.mnew = f.mold
    while (!converged) {
      
      #Prepare to iterate conjugate gradient descent:
      i=0
      f.old = log.lik(data, xx, event)
      converged.cg = FALSE
      while(!converged.cg) {
        i = i+1
        if(i == 1)
          t = 1
        else
          t = t/0.9/0.9
        
        dir.new = score(data, xx, event, specific = specific)
        
        dir.new = dir.new/sqrt(sum(dir.new^2))
        
        #Use conjugate gradient here
        if (i > 1) {
          cdir.old = cdir.old*max(0, sum(dir.new*(dir.new-dir.old))/(sum(dir.old^2)))
          dir.old = dir.new
          dir.new = dir.new + cdir.old
        } else {
          dir.old = dir.new
        }
        cdir.old = dir.new
        
        newxx = xx + dir.new * t
        
        foundStep = FALSE
        f.proposed = log.lik(data, newxx, event)
        
        #First make sure that we are at least improving
        while(!foundStep && !converged.cg){
          if(f.proposed > f.old)
            foundStep = TRUE
          else {
            t <- t*0.9
            if(t <= .Machine$double.eps){
              #Step size too small, just take original parameters
              converged.cg = TRUE
              newxx = xx
              cat(paste("Step size too small, converged, t = ", t, "\n")) 
              # cat(paste("Old likelihood: ", f.old, "new likelihood: ", f.proposed, "\n"))
            } else {
              newxx <- xx + dir.new * t
              f.proposed = log.lik(data, newxx, event)
            }
          }
        }
        #Now we try to find the optimal step
        foundStep = FALSE
        f.best = f.proposed
        while(!foundStep && !converged.cg){
          t <- t * 0.9
          if(t <= .Machine$double.eps){
            converged.cg = TRUE
            cat(paste("Step size too small, converged, t = ", t, "\n")) 
          }
          newxx <- xx + dir.new * t
          
          f.proposed <- log.lik(data, newxx, event)
          
          if(f.proposed > f.best){
            f.best = f.proposed
          }
          else{
            t <- t/0.9
            newxx <- xx + dir.new * t
            foundStep <- TRUE
          }
        }
        
        xx = newxx
        f.proposed = log.lik(data, xx, event)
        
        if(i%%100 == 0){
          cat(paste(i," iterations of CG, log-likelihood: ", f.proposed, "\n"))
          if(i%%1000 == 0)
            print.table(xx)
        }
        
        if((f.proposed - f.old) < f.old * tol){
          converged.cg = TRUE
          cat(paste("CG Converged, log-likelihood: ", f.proposed, "\n"))
          print.table(xx)
        }
        f.old = f.proposed
      }
      
      converged.irls = FALSE
      i = 0
      while(!converged.irls){
        #Do IRLS now
        i = i+1
        gammasigma = irls(data, xx, event)
        newxx = xx
        newxx[(p*d+p+1):(p*d+p+n)] = gammasigma
        # print.table(xx)
        
        f.proposed = log.lik(data, newxx, event)
        if(is.finite(f.proposed)){
          if((f.proposed - f.old) <= 0){
            converged.irls = TRUE
          }
        }
        xx = newxx
        
        cat(paste("log-likelihood after ", i, " iterations of IRLS: ", f.proposed, "\n"))
        f.old = f.proposed
      }
      
      f.mnew = f.proposed
      if ((f.mnew - f.mold) < f.mold * tol){
        converged = TRUE
        cat("M-step complete \n")
      }
      cat(paste("log-likelihood change: ", f.mnew - f.mold, "\n"))
      cat(paste("tolerance: ", f.mold * tol, "\n"))
      f.mold = f.mnew
    }
    
    d = length(unique(event))
    p = ncol(data)
    
    alpha = xx[1:(p*d)]
    beta = xx[(p*d+1):(p*d+p)]
    gamma = xx[(p*d+1+p):(p*d+p+n)]
    sigma = xx[p*d+p+n+1]
    data.new = E.step(alpha, beta, gamma, data, min.detect, event)
    
    check.old = check
    check = sum((data.new - data)^2, na.rm=TRUE) / sum(data^2, na.rm=TRUE)
    data = data.new
    
    cat(paste("E-step check: ", abs(check.old-check), "compared to ", tol, "\n"))
    if (abs(check.old - check) < tol) {
      if (tol<=sqrt(.Machine$double.eps))
        finished = TRUE
      tol = max(tol/2, sqrt(.Machine$double.eps))
      cat(paste("Iterating with tol=", tol, "\n", sep=""))
    }
  }
  
  #Compile the results and return
  result = list()
  result$data = data
  result$min.detect = min.detect
  result$event = event
  result$specific = specific
  result$alpha = matrix(alpha, d, p, byrow=TRUE)
  colnames(result$alpha) = names(result$beta)
  result$beta = beta
  result$gamma = gamma
  result$sigma = sigma
  
  return(result)
}