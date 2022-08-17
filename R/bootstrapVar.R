#' A boostap method to calculate the standard error of the variance
#'
#' This function uses the bootrstap method to calculate the uncertainty of the variance of a
#' given sample based on random resampling. The number of the resamples is a parameter (default
#' is 1000). Given that the resampling methods underestimate the uncertainty and thus
#' provide a biased estimation, we offer the the unbiased method as a default, although the user
#' may change this option through the biased parameter for experimental purposes (they are
#' strongly advised not to do for real problems with small samples).
#' @param biased A logical parameter to indicate if the user wants the biased version. Resampling
#' techniques always underestimate statistics like the variance or the standard error of it
#'  for small samples.
#' @param nrep The number of bootstrap resamples. Default is 1000. The higher the number of
#' the samples, the better the bootstrap outcome.
#' @inheritParams analyticVar
#' @return The standard error of the variance of \code{data}.
#' @keywords bootstrap uncertainty heteroplasmy resampling
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
#' bootstrapVar(mouseData1)
#'
#' # use the package data and load it to variable mouseData
#' mouseData=mousedataLE
#' # calculate the standard error of the variance for the LE oocyte sample #3
#' bootstrapVar(mouseData[,3])

bootstrapVar <- function(data,nrep=1000,biased=FALSE){
  if (typeof(data)!="double" || typeof(biased)!="logical" || (typeof(nrep)!="integer" && typeof(nrep)!="double")) {
    stop("Invalid data type(s). Check if the arguments' types are correct.")
  }
  if (length(data[which(is.na(data[]))])>length(data[which(!is.na(data[]))])) {
    warning("There were NA values in the input data and they were ommitted.")
  }
  if (nrep==0) {
    warning("The number of bootstrap cannot be 0: Changed to 1000.")
    nrep=1000
  }
  X.1 = data
  X.1=X.1[which(!is.na(X.1[]))]
  n=length(X.1)
  if (n<1) {
    stop("Invalid length of input data. It should be >0.")
  }
  #X.1=X.1/100
  h0=mean(X.1)

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

  return(heteropVar)
}
