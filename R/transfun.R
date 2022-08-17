#' A Function to cast real numbers to \code{[0,1]}
#'
#' A transformation function to cast any real number onto the interval \code{[0,1]}.
#' Equivalent to the inverse logit transform.
#' @param x a real number to be transformed.
#' @return The transformed cast of the input value to the interval \code{[0,1]}.
#' @keywords reverse logit
#' @export
#' @examples
#' transfun(5.1)

transfun = function(x) {
  return(1/(1+exp(-x)))
}
