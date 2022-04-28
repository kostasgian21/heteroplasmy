#' Calculate negative log likelihood for heteroplasmy measurements
#'
#' calculate negative log likelihood for a given set of heteroplasmy measurements
#' h and parameters theta = {logit(p), logit(b)} (we write h0 for p) we can do this
#' enforcing a particular h0 value (passed as an argument) or treating h0 as a fit
#' parameter (default) the logit transform is used to ensure h0 and b remain in
#' the \code{[0,1]} interval regardless of what real-valued argument the numerical optimiser
#' attempts
#' @param h0 Logical parameter. A particular h0 value  Default is to treat h0 as a fit parameter
#' @param h The heteroplasmy measurements.
#' @param theta Kimura parameters p (or h0 here) and b.
#' @return The negative log likelihood  for the input.
#' @keywords negative log likelihood kimura
#' @export
#' @examples
#'  X.1 = rnorm(50,0.5,0.1)
#' kimura_neg_loglik(c(0.5,0.91),X.1)

kimura_neg_loglik = function(theta, h, h0=F) {
  # get kimura parameters from argument
  b = transfun(theta[1])
  # if we haven't provided a specific h0, retrieve it as a parameter (only difference is the transformation step)
  if(h0 == F) {  h0 = transfun(theta[2]) }
  return(-sum(log(kimura::dkimura(h, h0, b))))
}
