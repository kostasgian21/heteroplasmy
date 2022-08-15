#' The inverse function to transfun
#'
#' A transformation function to inveret the effect of transfun (ie, the inverse logit transform).
#' @param x an integer number to be transformed.
#' @return The inveresed treansformation of transfun.
#' @keywords inverse logit
#' @export
#' @examples
#' invtransfun(0.71)

invtransfun = function(x) {
  return(log(x / (1-x)))
}
