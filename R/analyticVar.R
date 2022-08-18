#' Analytic calculation of the standard error of the variance
#'
#' This model-free method to calculate analytically the standard error of the variance.
#' It's based on the use of the appropriate h-statistic as an estimator, as default. It is offers a corrected version of the method described in Wonnapinij et al.. If instead the Wonnapinij et al. mehod shall be used, make method="Wonnapinij".
#' @param method What method to use for the estimation of the standard error of the variance. Accepted values are "hstatistic" (default) and "Wonnapinij".
#' @param normal Parameter that indicates if the normal approximation should be used
#' instead of the general formula from (Wilks, S. S. (1962).Mathematical Statistics).
#' Default is FALSE.
#' @param data The input data in the form of a vector. NA values are omitted.
#' @return The analytically derived standard error of the variance of \code{data}.
#' @keywords standard error variance heteroplasmy h-statistic
#' @export
#' @examples
#' # size of the sample
#' n=50
#' #generate a random sample of size n from a normal distribution
#' data_ex=rnorm(n,0.5,0.1)
#' analyticVar(data)
#'
#' mouseData=readHeteroplasmyData("HB")
#' mouseData1 = mouseData[which(!is.na(mouseData[,1])),1]
#' analyticVar(mouseData1,method="hstatistic")
#'
#' # use the package data and load it to variable mouseData
#' mouseData=mousedataLE
#' # calculate the standard error of the variance for the LE oocyte sample #3
#' bootstrapVar(mouseData[,3])



analyticVar <- function(data,normal=FALSE,method="hstatistic") {
  if (typeof(data)!="double" || typeof(normal)!="logical") {
    stop("Invalid data type(s). Check if the arguments' types are correct.")
  }
  if (method!="hstatistic" && method!="Wonnapinij") {
    warning("Invalid method. Please choose either 'hstatistic' (default value if empty) or 'Wonnapinij'. Otherwise default is used.")
  }
  if (length(data[which(is.na(data[]))])>length(data[which(!is.na(data[]))])) {
    warning("There were NA values in the input data and they were ommitted.")
  }
  X.1 = data
  X.1=X.1[which(!is.na(X.1[]))]
  X.1=as.vector(X.1)
  n=length(X.1)
  if (n<1) {
    stop("Invalid length of input data. It should be >0.0")
  }
  #X.1=as.data.frame(X.1)

  #X.1=X.1/100
  h0=mean(X.1)
  sampleVar=var(X.1)

  # SEW1 and SEW are the standard error of the variance SEW1 is the normal approximation and SEW2
  # is the generic type
  mu2=(1/n)*sum((X.1-h0)^2)
  mu4=(1/n)*sum((X.1-h0)^4)
  SEW1=sampleVar*sqrt(2/(n-1))# standard error,normal approximation

  if (method=="hstatistic") {
    D4=-3 * mu2^2 * (2 * n - 3) * n/((n - 1) * (n - 2) * (n - 3)) +
      (n^2 - 2 * n + 3) * mu4 * n/((n - 1) * (n - 2) * (n - 3))
  }
  if (method=="Wonnapinij") {
    D4=-(n-1)/n^3 * ((n^2-3*n+3)*mu4 + 3*(2*n-3)*mu2^2)
  }

  SEW2=sqrt(abs((1/n)*(D4-sampleVar^2*((n-3)/(n-1)))))
  if (normal==TRUE) {
    message("Normal approximation of the standard error of the variance")
    return(SEW1)
  }else{
    message(paste("Standard error of the variance using the",method,"method"))
    return(SEW2)
  }
}
