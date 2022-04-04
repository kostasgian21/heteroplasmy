#' Transformed heteroplasmy shift
#'
#' A numerical transformation of the heteroplasmy samples in order to work with the
#' heteroplasmy shifts across diverse samples (e.g., due to time or different tissue samples). This
#' transformation is used for comparing a heteroplasmy observation \emph{h} to a reference value
#' \ifelse{html}{\emph{h}\out{<sub>0</sub>}}{\eqn{h_0}}.
#' It corresponds to the formula:\cr \cr
#' \ifelse{html}{\eqn{\Delta}\out{h =ln((h(h<sub>0</sub> - 1))/(h<sub>0</sub>(h - 1)))}}{\deqn{\Delta h =  \ln \left( \frac{h (h_0 - 1)}{h_0 (h - 1)}\right)}}
#' @param h0 The reference heteroplasmy value. Should be in \code{[}0,1\code{/]}.
#' @param h The heteroplasmy observation. Can be either a single value or a vector of observations.
#' Every observation should be in \[0,1\] and [[0,1]].
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
#' heteroplasmyShift(data)
#'
#' mouseData=readHeteroplasmyData("HB")
#' mouseData1 = mouseData[which(!is.na(mouseData[,1])),1]
#' heteroplasmyShift(mouseData1,nrep=10000)

heteroplasmyShift <- function(h,h0) {
  if (typeof(h)!="double" || typeof(h0)!="double") {
    stop("Invalid data type(s). Check if the arguments' types are correct.")
  }
  if (length(h)==1) {
    if (h<0 || h0<0 || h>1 || h0>1) {
      stop("h or h0 should be in [0,1]")
    }
    deltaH=log((h*(h0-1))/(h0*(h-1)))
    warning("Input h type: single value")
    return(deltaH)
  }else{
    deltaH=log((h*(h0-1))/(h0*(h-1)))
    warning("Input h type: vector")
    return(deltaH)
  }
}
