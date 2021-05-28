#' A boostap method to calculate the standard error of the variance
#'
#' This function uses the bootrstap method to calculate the uncertainty of the variance of a
#' given sample based on random resampling. The number of the resamples is a parameter (default
#' is 1000) and along with the the "vanilla" version, we offer an optimized variation (using
#' the sigmaOpt parameter) which has been seen to improve the precision of the calculation (see
#' our report/paper). Given that the resampling methods underestimate the uncertainty and thus
#' provide a biased estimation, we offer the the unbiased method as a default, although the user
#' may change this option through the biased parameter for experimental purposes (they are
#' strongly advised not to do for real problems with small samples).
#' @param sigmaOpt The outcome of the bootstrap resampling with a fitted sigmoid function g(x)
#' with four parameters. Derived through simulations on both real heteroplasmy data and various
#' synthetic ones. Try the plotStdErrVar function in this package to observe it.
#' @param corrected Simle correction with a factor of 2.61 that was experimentally found. It is
#' included also in the case of sigmaOpt=TRUE
#' @param biased A logical parameter to indicate if the user wants the biased version. Resampling
#' techniques always underestimate statistics like the variance or the standard error of it
#'  for small samples.
#' @param nrep The number of bootstrap resamples. Default is 1000. The higher the number of
#' the samples, the better the bootstrap outcome.
#' @param data The input data in the form of a dataframe or matrix (which will be transformed into
#' a dataframe).
#' @keywords bootstrap,uncertainty,heteroplasmy,resampling
#' @export
#' @examples
#' # size of the sample
#' n=50
#' #generate a random sample of size n from a normal distribution
#' data_ex=rnorm(n,0.5,0.1)
#' bootstrapVar(data)
#'
#' mouseData=readHeteroplasmyData("HB")
#' mouseData1 = mouseData[which(!is.na(mouseData[,1])),1]
#' bootstrapVar(mouseData1,sigmaOpt=TRUE)

bootstrapVar <- function(data,nrep=1000,biased=FALSE,corrected=FALSE,sigmaOpt=FALSE) {
  n=length(which(!is.na(data)))
  X.1 = data
  #X.1=X.1/100
  h0=mean(X.1)
  opti = 2.61
  coeB=-4.1434
  coeC=-0.7387
  coeD=1.0125
  coeE=2.8390

  X.1_boot=as.data.frame(X.1)
  boots = vector()
  bootVars=vector()
  bootVarsBiased=vector()

  for (i in 1:nrep) {
    repeat {
      boot= sample(1:nrow(X.1_boot), nrow(X.1_boot), replace=TRUE)
      boots=X.1_boot[boot,1]
      if (!is.na(var(boots,na.rm = TRUE))) {
        break
      }
    }
    #append the corrected variance of the i-th bootstrap sample
    bootVars= c(bootVars,var(boots,na.rm = TRUE)*(n/(n - 1)))
    bootVarsBiased= c(bootVarsBiased,var(boots,na.rm = TRUE))

  }

  # standard error of the bootstrap samples
  SEB=sd(bootVars)
  heteropVar=SEB
  if (biased==TRUE) {
    warning("Cautious, biased calculation always underestimates the variance!")
    SEB=sd(bootVarsBiased)
    heteropVar=SEB
  }

  if (corrected==TRUE) {
    message("Correceted by a constant factor")
    heteropVar=SEB*((n - opti)/n)
  }
  if (sigmaOpt==TRUE) {
    message("Correceted by a fitted sigmoid function, ie a four-parameter log-logistic function")
    # heteropVar=SEB*((n - opti)/n)/(coeC+(coeD-coeC)/(1+exp(coeB*(log(n)-log(coeE)))))
    heteropVar=SEB*(((n - opti)*(1+exp(coeB*(log(n)-log(coeE)))))/(n*(coeC+(coeD-coeC))))

  }

  return(heteropVar)
}
