#' Max likelihood estimation for Kimura distribution parameters
#'
#' Using maximum likelihood to estimate the parameters of a fitted Kimura distribution to the input sample values.
#' @inheritParams maxlik
#' @return The maximum likelihood estimates for a fitted Kimura distribution parameters.
#' @author Kostas and Iain, \email{us@@example.com}
#' @references \href{http://example.com}{Site or paper}
#' @seealso \code{\link{readHeteroplasmyData}}
#' @keywords heteroplasmy,maximum likelihood
#' @export
#' @examples
#' # size of the sample
#' n=50
#' #generate a random sample of size n from a normal distribution
#' data_ex=rnorm(n,0.5,0.1)
#' estimate_parameters_ml(data_ex)
#'
#' mouseData=readHeteroplasmyData("LE")
#' mouseData1 = mouseData[which(!is.na(mouseData[,1])),1]
#' estimate_parameters_ml(mouseData1)

estimate_parameters_ml = function(h) {
  # find best transformed parameters h0 and b
  best = stats::optim(c(0.5, 0.5), kimura_neg_loglik, h=h, h0=F, hessian=F)

  best$b.hat = transfun(best$par[1])
  best$h0.hat = transfun(best$par[2])

  return(c(best$h0.hat, best$b.hat))
}
