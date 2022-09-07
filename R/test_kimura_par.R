#' Generalised MC KS test for genetic drift
#'
#' This function is a generalisation of the test_kimura function from the lbozhilova/kimura package. It corresponds to a Monte Carlo Kolmogorov-Smirnov test to detect genetic drift by examining deviation from a Kimura distribution. (add references)
#'
#' @inheritParams estimate_parameters_ml
#' @param p The p parameter of the Kimura distribution. Should be in \code{[0,1]}.
#' @param b The b parameter of of the Kimura distribution.  Should be in \code{[0,1]}.
#' @param num_MC number of Monte Carlo runs
#' @param round a logical argument. True if heteroplasmy fractions are rounded to two significant digits.
#'
#' @return object of class htest
#' @export
#'
#' @examples
#' data_ex=rnorm(n,0.5,0.1)
#'   fit = estimate_parameters_ml(data_ex)
#'   p=fit[1]
#'   b=fit[2]
#' test_kimura_par(data_ex,p,b)
test_kimura_par <- function(h, p, b, num_MC = 1000, round = TRUE) {
  cdf_kimura <- kimura::.pkimura_full(p, b)
  h_mat <- cbind(h, matrix(kimura::rkimura(num_MC * length(h), p, b), ncol = num_MC))
  if (round)
    h_mat <- round(h_mat, 2)
  get_ks_statistic <- function(h) {
    ecdf_h <- stats::ecdf(h)(seq(0, 1, 1e-4))
    max(abs(ecdf_h - cdf_kimura))
  }
  D <- unlist(apply(h_mat, 2, get_ks_statistic))
  D <- unname(D)
  output <- list(
    statistic = c("D" = D[1], "p" = p, "b" = b),
    p.value = unname(sum(D >= D[1]) / length(D)),
    alternative = "one-sided",
    method = "Monte Carlo Kolmogorov-Smirnov",
    data.name = paste0(deparse(substitute(h)),
                       " and Kimura(", round(p, 4), ", ", round(b, 4), ")"),
    D_sample = D)
  class(output) <- "htest"
  output
}
