#' Joint negative log likelihood function for several heteroplasmy measurements
#'
#' joint negative log likelihood function for several families' heteroplasmy
#' measurements
#' theta = \code{[ b, h0.1, h0.2, ... ]} (use h values if use.h0s=F, otherwise initial heteroplasmies are enforced via h0s).
#' @param h0s TBDD
#' @param use.h0s Logical parameter. TBD
#' @param hlist TBDD list of different sets of heteroplasmy measurements
#' @inheritParams kimura_neg_loglik
#' @return The negative log likelihood  for the list of inputs.
#' @keywords joint negative log likelihood kimura
#' @export
#' @examples
#'  X.1 = rnorm(50,0.5,0.1)
#' joint_neg_log_lik(c(0.5,0.91),X.1)

joint_neg_log_lik = function(theta, hlist, use.h0s=F, h0s=-1) {
  llik = 0
  for(i in 1:length(hlist)) {
    if(use.h0s == F) {
      llik = llik + kimura_neg_loglik(c(theta[1], theta[i+1]), hlist[[i]], h0=F)
    } else {
      llik = llik + kimura_neg_loglik(theta[1], hlist[[i]], h0=h0s[i])
    }
  }
  return(llik)
}
