#' Bootstrap estimates for parameters and confidence intervals for heteroplasmy data
#'
#' compute bootstrap estimates for parameters and confidence intervals for a given heteroplasmy set
#' Wwe can do this while imposing a specific h0 as an argument or allowing a search over h0 values
#' @param nboot The number of bootstrap samples. Default value is 1000
#' @inheritParams maxlik
#' @return The maximum likelihood for the input data according to the Kimura distribution (using bootstrapping)
#' @keywords joint negative log likelihood kimura
#' @export
#' @examples
#'  X.1 = rnorm(50,0.5,0.1)
#' maxlik(X.1,nboot=10000,0.95)

maxlikboot = function(h, nboot=1000, conf.level = 0.95, h0=F) {
  # if we have enforced a particular h0
  if(h0 != F) {
    boot.b = c()
    # loop over bootstrap resamples
    for(boot in 1:nboot) {
      # construct bootstrap sample
      hboot = sample(h, replace=T)
      # find best transformed parameter b
      boot.best = optim(c(0.5), kimura_neg_loglik, h=hboot, h0=h0, method="Brent", lower=-30, upper=30)
      # record back-transformed parameter for this resample
      boot.b = c(boot.b, transfun(boot.best$par[1]))
    }
    # do the optimisation for the non-resampled set
    best = optim(c(0.5), kimura_neg_loglik, h=h, h0=h0)
    best$b.hat = transfun(best$par[1])
    best$n.hat = 1/(1-best$b.hat)
    best$h0.hat = h0

    # get stats and confidence intervals from the bootstrap distribution
    best$b.bhat = mean(boot.b)
    best$h0.ci = c(h0,h0)
    conf.1 = (1-conf.level)/2
    conf.2 = 1-conf.1
    best$b.ci = c(quantile(boot.b, conf.1), quantile(boot.b, conf.2))

    best$n.bhat = 1/(1-best$b.bhat)
    best$n.ci = 1/(1-best$b.ci)
  } else {
    boot.h0 = boot.b = c()
    # loop over bootstrap resamples
    for(boot in 1:nboot) {
      hboot = sample(h, replace=T)
      # find best transformed parameters h0 and b
      boot.best = optim(c(0.5, 0.5), kimura_neg_loglik, h=hboot, h0=F)
      # record back-transformed parameters for this resample
      boot.h0 = c(boot.h0, transfun(boot.best$par[2]))
      boot.b = c(boot.b, transfun(boot.best$par[1]))
    }
    # do the optimisation for the non-resampled set
    best = optim(c(0.5, 0.5), kimura_neg_loglik, h=h)
    best$h0.hat = transfun(best$par[2])
    best$b.hat = transfun(best$par[1])
    best$n.hat = 1/(1-best$b.hat)

    # get stats and confidence intervals from the bootstrap distribution
    best$h0.bhat = mean(boot.h0)
    best$b.bhat = mean(boot.b)
    best$h0.ci = c(quantile(boot.h0, 1-conf.level), quantile(boot.h0, conf.level))
    best$b.ci = c(quantile(boot.b, 1-conf.level), quantile(boot.b, conf.level))

    best$n.bhat = 1/(1-best$b.bhat)
    best$n.ci = 1/(1-best$b.ci)
  }
  return(best)
}
