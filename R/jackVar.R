#' A jackknife method to compute the uncertainty of  heteroplasmy data
#'
#' Similarly to the main \emph{bootstrapVar} function that implements the bootstrap method to measure the
#' standard error of the variance, the jackknife technique is another resampling method that can
#' be used for the same purpose. Unlike \emph{bootstrapVar}, \emph{jackVar} (and every jackknife-like method) is deterministic and does not rely on randomness, but instead it uses removals of the sample points, one each
#' time to calculate different sub-samples of size (n-1). Note that the size of the input data
#' should be strictly greater than 1.
#' @inheritParams analyticVar
#' @return The analytically derived standard error of the variance of \code{data}  and the mean of the bootsrap samples means..
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
#'
#' # use the package data and load it to variable mouseData
#' mouseData=mousedataLE
#' # calculate the standard error of the variance for the LE oocyte sample #3
#' bootstrapVar(mouseData[,3])
#' \dontrun{
#' #input data of size 1 will fail
#' data_ex=rnorm(1,0.5,0.1)
#' jackVar(data)
#' }

jackVar <- function(data) {
  if (typeof(data)!="double") {
    stop("Invalid data type. Check if the arguments' types are correct.")
  }
  if (length(data[which(is.na(data[]))])>length(data[which(!is.na(data[]))])) {
    warning("There were NA values in the input data and they were ommitted.")
  }
  X.1 = data
  X.1=X.1[which(!is.na(X.1[]))]
  n=length(X.1)
  if (n<2) {
    stop("Invalid length of input data for jackknife. By definition, it should be >1.")
  }
  #X.1=X.1/100
  h0=mean(X.1)
  if (n<2) {
    stop("Input data should have size >1")
    return()
  }

  # jackknife resampling -- pseudo reports the vars of each subsample
  # and pseudoSE is for the calculation of its standard error
  pseudoV <- numeric(length(X.1))
  pseudoM <- numeric(length(X.1))
  for (ww in 1:length(X.1)) {
    pseudoV[ww] <- var(X.1[-ww])
    pseudoM[ww] <- mean(X.1[-ww])
  }
  pseudoSE <- sqrt(((n - 1)/n) * sum((pseudoV - mean(pseudoV))^2))
  meanJack=mean(pseudoM)
  return(c(pseudoSE, meanJack))
}
