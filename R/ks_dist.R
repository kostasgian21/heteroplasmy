#' Kolmogorov-Smirnov distance function
#'
#' A function to calculate the Kolmogorov-Smirnov distance. It is used in estimate_parameters_ks to estimate the parameterr that minimize the distance in the optim function.
#' @param ecdf A vector containing 10000 values from the empirical cumulative distribution function of the input heteroplasmy data vector. To be used in the optim function in estimate_parameters_ks.
#' @param theta A vector. with two elements. The two Kimura parameters h0 and b.
#' @return The maximum distance from the theoretical Kimura distribution.
#' @keywords Kolmogorov Smirnov distance
#' @export
#' @examples
#' h=rnorm(20,0.5,0.1)
#' ecdf_h <- (stats::ecdf(h))(seq(0, 1, 1e-04))
#' ks_dist(c(0.5,0.95),ecdf_h)

ks_dist = function(theta, ecdf) {
  b = transfun(theta[2])
  h0 = transfun(theta[1])
  cdf_kimura <- kimura::.pkimura_full(h0, b)
  return(max(abs(ecdf - cdf_kimura)))
}
