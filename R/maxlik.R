#' compute maximum likelihood parameters and confidence intervals for heteroplasmy data
#'
#' compute maximum likelihood parameters and confidence intervals for a given heteroplasmy set.
#' We can do this while imposing a specific h0 as an argument or allowing a search over h0 values.
#' @param h0 Logical parameter. A particular h0 value  Default is to treat h0 as a fit parameter
#' @param conf.level The preferred confidence interval calculation, Default value is 0,95 (95%).
#' @param h The heteroplasmy measurements.
#' @return The maximum likelihood for the input data according to the Kimura distribution
#' @keywords joint negative log likelihood kimura
#' @export
#' @examples
#'  X.1 = rnorm(50,0.5,0.1)
#' maxlik(X.1,0.95)

maxlik = function(h, conf.level = 0.95, h0 = F) {
  # if we have enforced a particular h0
  if(h0 != F) {
    # find best transformed parameter b
    best = stats::optim(c(0.5), kimura_neg_loglik, h=h, h0=h0, method="Brent", lower=-30, upper=30, hessian=T)
    # add back-transformed parameters to return structure
    # h0 is fixed here; n is straightforward function of b
    best$h0.hat = h0
    best$b.hat = transfun(best$par[1])
    best$n.hat = 1/(1-best$b.hat)

    # the hessian from the optimisation process is the fisher information matrix.
    # its inverse (via "solve") can be used to construct confidence intervals on our parameter(s)
    fisher.matrix = best$hessian
    crit = qnorm((1 + conf.level)/2)
    inv.fisher.matrix <- solve(fisher.matrix)
    b.ci = best$par[1] + c(-1, 1) * crit * sqrt(inv.fisher.matrix[1, 1])
    best$h0.ci = c(h0, h0)
    best$b.ci = transfun(b.ci)

    best$n.hat = 1/(1-best$b.hat)
    best$n.ci = 1/(1-best$b.ci)
  } else {
    # find best transformed parameters h0 and b
    best = optim(c(0.5, 0.5), kimura_neg_loglik, h=h, h0=F, hessian=T)
    # add back-transformed parameters to return structure
    best$b.hat = transfun(best$par[1])
    best$h0.hat = transfun(best$par[2])
    best$n.hat = 1/(1-best$b.hat)

    # the hessian from the optimisation process is the fisher information matrix.
    # its inverse (via "solve") can be used to construct confidence intervals on our parameter(s)
    fisher.matrix = best$hessian
    crit = qnorm((1 + conf.level)/2)
    inv.fisher.matrix <- solve(fisher.matrix)
    b.ci = best$par[1] + c(-1, 1) * crit * sqrt(inv.fisher.matrix[1, 1])
    h0.ci = best$par[2] + c(-1, 1) * crit * sqrt(inv.fisher.matrix[2, 2])
    best$b.ci = transfun(b.ci)
    best$h0.ci = transfun(h0.ci)

    best$n.hat = 1/(1-best$b.hat)
    best$n.ci = 1/(1-best$b.ci)
  }

  return(best)
}
