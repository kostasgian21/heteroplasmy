#' The bootstrapVar function with default correction
#'
#' This function is simply the bootstrapVar with the correction argument being TRUE. it is
#' provided as a seperate function for usability.The function can be used beyond heteroplasmy data,
#' therefore one can use it to calculate the standard error of the variance for samples where other
#' approaches may not fit, eg when the sample size is too small and/or the population distribution
#' is not Gaussian (or not known at all).
#' @param nrep The number of bootstrap resamples. Default is 1000. The higher the number of
#' the samples, the better the bootstrap outcome (see \code{\link[graphics]{par}}).
#' @param data The input data in the form of a dataframe or matrix (which will be transformed into
#' a dataframe). NA values are omitted.
#' @return The analytically derived standard error of the variance of \code{data}.
#' @author Kostas and Iain, \email{us@@example.com}
#' @references \url{https://en.wikipedia.org/}
#' @seealso \code{\link{bootstrapVar}}
#' @keywords bootstrapVarCor, fitted
#' @export
#' @examples
#' # size of the sample
#' n=50
#' #generate a random sample of size n from a normal distribution
#' data_ex=rnorm(n,0.5,0.1)
#' bootstrapVarCor(data)
#'
#' mouseData=readHeteroplasmyData("HB")
#' mouseData1 = mouseData[which(!is.na(mouseData[,1])),1]
#' bootstrapVarCor(mouseData1,nrep=10000)

bootstrapVarCor <- function(data,nrep=1000) {
  rep=nrep
  heteropVar=bootstrapVar(data,nrep=rep)
  return(heteropVar)
}
