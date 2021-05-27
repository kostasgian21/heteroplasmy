#' A Function to read mouse heteroplasmy data
#'
#' This function allows you to read mouse heteroplasmy data from external files.
#' @param normal Parameter that indicates if the normal approximation should be used
#' instead of the general formula from (Wilks, S. S. (1962).Mathematical Statistics).
#' Default is FALSE.
#' @param data The input data in the form of a dataframe or matrix (which be transformed into
#' a dataframe).
#' @keywords standard error of the variance, heteroplasmy
#' @export
#' @examples
#' analyticVar(data)


analyticVar <- function(data,normal=FALSE) {
  n=length(which(!is.na(data)))
  X.1=as.data.frame(X.1)
  X.1 = data
  #X.1=X.1/100
  h0=mean(X.1)
  sampleVar=var(X.1)

  # SEW1 and SEW are the standard error of the variance from
  # Wonnapinij et al. SEW1 is the normal approximation and SEW2
  # is the generic type
  mu2=(1/n)*sum((X.1-h0)^2)
  mu4=(1/n)*sum((X.1-h0)^4)
  D4=((n-1)/(n^3))*((n^2-3*n+3)*mu4+3*(2*n-3)*mu2^2)
  SEW1=sampleVar*sqrt(2/(n-1))# standard error,normal approximation
  SEW2=sqrt((1/n)*(D4-sampleVar^2*((n-3)/(n-1))))# standard error, not normal structure
  if (normal==TRUE) {
    message("Normal approximation of the standard error of the variance")
    return(SEW1)
  }else{
    return(SEW2)
  }
}
