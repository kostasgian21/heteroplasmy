#' Max likelihood estimation of Kimura parameters that minimize the KS statistic
#'
#' Using maximum likelihood to estimate the parameters of a fitted Kimura distribution to the input sample values that minimizes the KS statistic. It is used to showcase that using the KS statistic to support selection needs extra caution.
#' @inheritParams estimate_parameters_ml
#' @return The maximum likelihood estimates for a fitted Kimura distribution parameters that minimze the KS statistic of a KS test.
#' @author Kostas and Iain, \email{us@@example.com}
#' @references \href{http://example.com}{Site or paper}
#' @seealso \code{\link{readHeteroplasmyData}}
#' @keywords heteroplasmy,maximum likelihood
#' @export
#’ @importFrom kimura
#' @examples
#' # size of the sample
#' n=50
#' #generate a random sample of size n from a normal distribution
#' data_ex=rnorm(n,0.5,0.1)
#' estimate_parameters_ml(data_ex)
#'
#' mouseData=readHeteroplasmyData("LE")
#' mouseData1 = mouseData[which(!is.na(mouseData[,1])),1]
#' estimate_parameters_ks(mouseData1)

estimate_parameters_ks = function(h) {
  ecdf_h <- (stats::ecdf(h))(seq(0, 1, 1e-04))

  mom.best = kimura::estimate_parameters(h)
  mom.p = mom.best[1]
  mom.b = mom.best[2]
  ml.best = estimate_parameters_ml(h)
  ml.p = ml.best[1]
  ml.b = ml.best[2]

  mom.best = stats::optim(c(invtransfun(mom.p), invtransfun(mom.b)), ks_dist, ecdf = ecdf_h)
  ml.best = stats::optim(c(invtransfun(ml.p), invtransfun(ml.b)), ks_dist, ecdf = ecdf_h)

  if(ml.best$value < mom.best$value) {
    b.hat = transfun(ml.best$par[2])
    h0.hat = transfun(ml.best$par[1])
  } else {
    b.hat = transfun(mom.best$par[2])
    h0.hat = transfun(mom.best$par[1])
  }

  return(c(h0.hat, b.hat))
}
