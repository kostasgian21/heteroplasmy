#' Transformed heteroplasmy shift
#'
#' A numerical transformation of the heteroplasmy samples in order to work with the
#' heteroplasmy shifts across diverse samples (e.g., due to time or different tissue samples).
#' @param nrep The number of bootstrap resamples. Default is 1000. The higher the number of
#' the samples, the better the bootstrap outcome (see \code{\link[graphics]{par}}).
#' @param data The input data in the form of a dataframe or matrix (which will be transformed into
#' a dataframe).
#' @author Kostas and Iain, \email{us@@example.com}
#' @references \url{https://en.wikipedia.org/}
#' @seealso \code{\link{readHeteroplasmyData}}
#' @keywords heteroplasmy,transformation,shift
#' @export
#' @examples
#' # size of the sample
#' n=50
#' #generate a random sample of size n from a normal distribution
#' data_ex=rnorm(n,0.5,0.1)
#' heteroplsamyShift(data)
#'
#' mouseData=readHeteroplasmyData("HB")
#' mouseData1 = mouseData[which(!is.na(mouseData[,1])),1]
#' heteroplsamyShift(mouseData1,nrep=10000)

heteroplsamyShift <- function(h,h0) {
  deltaH=log((h*(h0-1))/(h0*(h-1)))
  return(deltaH)
}
