#' Likelihood ratio between two heteroplasmy samples
#'
#' A function to perform likelihood ratio test exploring difference in bottleneck size between two heteroplasmy samples.
#' @param h2.h0set Logical parameter.TBD
#' @param h1.h0set Logical parameter.TBD
#' @param use.h0s Logical parameter. TBD
#' @param h2 The two heteroplasmy sample.
#' @param h1 The first heteroplasmy sample.
#' @return The maximum likelihood for the input data according to the Kimura distribution (using bootstrapping)
#' @keywords summary statistic
#' @export
#' @examples
#'  X.1 = rnorm(50,0.5,0.1)
#'  X.2 = rnorm(50,0.5,0.1)
#' kimura_lrt(X.1,X.2)

kimura_lrt = function(h1, h2, use.h0s=F, h1.h0set=0, h2.h0set=0) {

  comp.df = data.frame()

  if(use.h0s) {
    # single-param optimisation for PT, MT, and combined h sets
    h1.best = stats::optim(-3, joint_neg_log_lik, hlist=h1, use.h0s=T, h0s=h1.h0set, method="Brent", lower=-30, upper=30, hessian=T)
    h2.best = stats::optim(-3, joint_neg_log_lik, hlist=h2, use.h0s=T, h0s=h2.h0set, method="Brent", lower=-30, upper=30, hessian=T)
    both.best = stats::optim(-3, joint_neg_log_lik, hlist=c(h1,h2), use.h0s=T, h0s=c(h1.h0set,h2.h0set), method="Brent", lower=-30, upper=30, hessian=T)
  } else {
    # multi-param optimisation for PT, MT, and combined h sets
    h1.best = stats::optim(c(-3, rep(0.5, length(h1))), joint_neg_log_lik, hlist=h1, use.h0s=F, hessian=T)
    h2.best = stats::optim(c(-3, rep(0.5, length(h2))), joint_neg_log_lik, hlist=h2, use.h0s=F, hessian=T)
    both.best = stats::optim(c(-3, rep(0.5, length(c(h1,h2)))), joint_neg_log_lik, hlist=c(h1, h2), use.h0s=F, hessian=T)
  }

  # use Fisher information to get confidence intervals on bottleneck size estimates
  conf.level = 0.95
  crit = qnorm((1 + conf.level)/2)
  h1.ci = h1.best$par[1] + c(-1, 1) * crit * sqrt(solve(h1.best$hessian)[1, 1])
  h2.ci = h2.best$par[1] + c(-1, 1) * crit * sqrt(solve(h2.best$hessian)[1, 1])
  both.ci = both.best$par[1] + c(-1, 1) * crit * sqrt(solve(both.best$hessian)[1, 1])

  h1.n.hat = 1/(1-transfun(h1.best$par[1]))
  h1.n.ci = 1/(1-transfun(h1.ci))
  h2.n.hat = 1/(1-transfun(h2.best$par[1]))
  h2.n.ci = 1/(1-transfun(h2.ci))
  both.n.hat = 1/(1-transfun(both.best$par[1]))
  both.n.ci = 1/(1-transfun(both.ci))

  # get log-likelihoods for the separate (MT =/= PT) and both (MT = PT) models (returned values are negative log likelihoods, so take negatives here)
  sep.llik = -h2.best$value-h1.best$value
  both.llik = -both.best$value

  # likelihood ratio test and p-value from chi-squared distribution (one dof -- one parameter difference)
  lrt = -2*(both.llik - sep.llik)
  pval = stats::pchisq(lrt, 1, lower.tail=F)
  comp.df = rbind(comp.df, data.frame(h1.n.hat, h1.n.ci[1], h1.n.ci[2], h2.n.hat, h2.n.ci[1], h2.n.ci[2], both.n.hat, both.n.ci[1], both.n.ci[2], pval))
  return(comp.df)
}
