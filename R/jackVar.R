#' A Function to read mouse heteroplasmy data
#'
#' This function allows you to read mouse heteroplasmy data from external files.
#' @param nameD Either "HB" or "LE".
#' @keywords heteroplasmy data
#' @export
#' @examples
#' jackVar(data)

jackVar <- function(data) {
  n=length(which(!is.na(data)))
  X.1 = data
  #X.1=X.1/100
  h0=mean(X.1)
  if (n<2) {
    stop("Input data should have size >1")
    return()
  }

  # jackknife resampling -- pseudo reports the vars of each subsample
  # and pseudoSE is for the calculation of its standard error
  pseudo <- numeric(length(X.1))
  for (ww in 1:length(X.1)) {
    pseudo[ww] <- var(X.1[-ww])
  }
  pseudoSE <- sqrt(((n - 1)/n) * sum((pseudo - mean(pseudo))^2))
  return(pseudoSE)
}
