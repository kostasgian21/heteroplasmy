#' The inverse function to transfun
#'
#' A transformation function to inveret the effect of transfun (ie, the inverse logit transform).
#' @param x a real number in \code{[0,1]} to be transformed into a real.
#' @return The inveresed treansformation of transfun.
#' @keywords inverse logit
#' @export
#' @examples
#' invtransfun(0.71)

invtransfun = function(x) {
  return(log(x / (1-x)))
}
