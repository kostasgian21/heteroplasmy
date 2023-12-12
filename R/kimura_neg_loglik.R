#' Calculate negative log-likelihood for heteroplasmy measurements
#'
#' Calculate negative log-likelihood for a given set of heteroplasmy measurements
#' \emph{h} and parameters theta = {logit(p), logit(b)} (we write \emph{h0} for \emph{p}). We can do this by
#' enforcing a particular \emph{h0} value (passed as an argument) or treating \emph{h0} as a fit
#' parameter (default). The logit transform is used to ensure \emph{h0} and \emph{b} remain in
#' the \code{[0,1]} interval regardless of what real-valued argument the numerical optimiser
#' attempts.
#' @inheritParams maxlik
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
