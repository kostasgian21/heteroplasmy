#' calculate various statistics for a heteroplasmy set h
#'
#' can enforce an initial h0 or leave as a free parameter. Can can use population or sample statistics. analyticVar offers a simplified version of this function to compute the standard error of the variance.
#' @param usepopn Logical parameter. Use of population or sample statistics (T and F, respectively)
#' @param h0 Logical parameter. A particular h0 value  Default is to treat h0 as a fit parameter
#' @param h The heteroplasmy measurements.
#' @return The maximum likelihood for the input data according to the Kimura distribution (using bootstrapping)
#' @keywords summary statistic
#' @export
#' @examples
#'  X.1 = rnorm(50,0.5,0.1)
#' hstats(X.1)

hstats = function(h, h0=F, usepopn=F) {
  # initialise list of results
  statres = list()
  # basic stats
  n = length(h)
  # if we're not enforcing a particular h0, compute summary stats
  if(h0 == F) {
    hbar = mean(h)
    s2 = var(h)
    if(usepopn == T) { s2 = (n-1)*s2/n }
    sehbar = sqrt(s2/n)
  } else {
    # otherwise, compute stats assuming the given h0
    hbar = h0
    if(usepopn == F) { s2 = 1/(n-1) * sum((h-hbar)**2) } else { s2 = 1/n * sum((h-hbar)**2) }
    sehbar = 0
  }
  # standard error on the variance and confidence intervals (Wonnapinij) - - cis rely on normality assumption of course
  sev = analyticVar(h)
  statres$moment.h0.hat = hbar
  statres$moment.h0.ci = c( hbar-1.96*sehbar, hbar+1.96*sehbar )

  # deal with the case where we have zero variance (identical samples). In practise this only happens in homoplasmic cases, where b and n are undefined (no information about bottleneck size)
  if(s2 == 0) {
    statres$moment.n.hat = NA
    statres$moment.n.ci = c(NA, NA)
  } else {
    statres$moment.n.hat = hbar*(1-hbar)/s2
    statres$moment.n.ci = c( (hbar*(1-hbar)/(s2+1.96*sev)), (hbar*(1-hbar)/(s2-1.96*sev)))
  }

  # for Kimura fit, decide whether to use bootstrapping or Fisher information for uncertainty
  # looks like homoplasmic samples reward bootstrapping (numerical convergence issues), and others prefer Fisher
  if(all(h == 0 | h == 1)) {
    print("Homoplasmic samples -- bootstrapping")
    fit = maxlikboot(h, h0=h0)
  } else {
    print("Heteroplasmic samples -- using Fisher information")
    fit = maxlik(h, h0=h0)
  }
  # estimates and confidence intervals from Kimura fit
  statres$fit.h0.hat = fit$h0.hat
  statres$fit.h0.ci = fit$h0.ci
  statres$fit.n.hat = fit$n.hat
  statres$fit.n.ci = fit$n.ci

  if(is.na(statres$moment.n.hat)) {
    statres$preferred = "fit"
  } else {
    # which approach is "better"? if confidence intervals on n go below 1, discard that approach; otherwise if one set of confidence intervals fits inside the other, favour that one
    if(statres$moment.n.ci[1] < 0 | statres$moment.n.ci[2] < 0 | (statres$fit.n.ci[1] > statres$moment.n.ci[1] & statres$fit.n.ci[2] < statres$moment.n.ci[2])) {
      statres$preferred = "fit"
    } else if(statres$fit.n.ci[1] < statres$moment.n.ci[1] & statres$fit.n.ci[2] > statres$moment.n.ci[2]) {
      statres$preferred = "moment"
    } else {
      statres$preferred = "neither"
    }
  }

  return(statres)
}
