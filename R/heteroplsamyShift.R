#' Transformed heteroplasmy shift
#'
#' A numerical transformation of the heteroplasmy samples in order to work with the
#' heteroplasmy shifts across diverse samples (e.g., due to time or different tissue samples). This
#' transformation is used for comparing a heteroplasmy observation \emph{h} to a reference value
#' \ifelse{html}{\emph{h}\out{<sub>0</sub>}}{\eqn{h_0}}.
#' It corresponds to the formula:\cr \cr
#' \ifelse{html}{\eqn{\Delta}\out{h =ln((h(h<sub>0</sub> - 1))/(h<sub>0</sub>(h - 1)))}}{\deqn{\Delta h =  \ln \left( \frac{h (h_0 - 1)}{h_0 (h - 1)}\right)}}
#' @param h0 The reference heteroplasmy value.
#' @param h The heteroplasmy observation.
#' @return The Transformed heteroplasmy shift.
#' @author Kostas and Iain, \email{us@@example.com}
#' @references \href{http://example.com}{Site or paper}
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
