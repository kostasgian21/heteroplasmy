#' A jackknife method to compute the uncertainty of  heteroplasmy data
#'
#' Similarly to the main bootstrapVar function that implements the bootstra method to measure the
#' standard error of the variance, the jackknife technique is another resampling method that can
#' be used for the same purpose. Unlike bootstrapVar, jackVar (and very jackknife method) is deterministic
#' and doesn not rely on randomness, but instead it uses removals of the sample points, one each
#' time to calculate different sub-samples of size (n-1). Note that the size of the input data
#' should be strictly greater than 1.
#' @param data The input data in the form of a dataframe or matrix (which will be transformed into
#' a dataframe). Its size should be >=2
#' @keywords jackknife uncertainty heteroplasmy resampling
#' @export
#' @examples
#' # size of the sample
#' n=50
#' #generate a random sample of size n from a normal distribution
#' data_ex=rnorm(n,0.5,0.1)
#' jackVar(data)
#'
#' mouseData=readHeteroplasmyData("HB")
#' mouseData1 = mouseData[which(!is.na(mouseData[,1])),1]
#' jackVar(mouseData1)

jackVar <- function(data) {
  n=length(which(!is.na(data)))
  X.1 = data
  #X.1=X.1/100
  h0=mean(X.1)
  if (n<2) {
    stop("Input data should have size >1")
    return()
  }

  # jackknife resampling -- pseudo reports the vars of each subsample
  # and pseudoSE is for the calculation of its standard error
  pseudo <- numeric(length(X.1))
  for (ww in 1:length(X.1)) {
    pseudo[ww] <- var(X.1[-ww])
  }
  pseudoSE <- sqrt(((n - 1)/n) * sum((pseudo - mean(pseudo))^2))
  return(pseudoSE)
}
