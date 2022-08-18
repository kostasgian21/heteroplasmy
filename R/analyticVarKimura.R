#' Kimura-based calculation of the standard error of the variance
#'
#' A Kimura-distribution-based method calculates analytically the standard error of the variance (see Eq. 11 in Wonnapinij et al., 2010).
#' It'uses only the Kimura distribution parameters p and b (passed as arguments).
#' @inheritParams analyticVar
#' @inheritParams test_kimura_par
#' @return The analytically derived standard error of the variance of \code{data}.
#' @keywords standard error variance heteroplasmy Kimura
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
#' analyticVar(mouseData1)
#'
#' # use the package data and load it to variable mouseData
#' mouseData=mousedataLE
#' # calculate the standard error of the variance for the LE oocyte sample #3
#' analyticVarKimura(mouseData[,3])



analyticVarKimura <- function(data,p,b) {
  if (typeof(data)!="double") {
    stop("Invalid data type(s). Check if the arguments' types are correct.")
  }
  if (length(data[which(is.na(data[]))])>length(data[which(!is.na(data[]))])) {
    warning("There were NA values in the input data and they were ommitted.")
  }
  X.1 = data
  X.1=X.1[which(!is.na(X.1[]))]
  X.1=as.vector(X.1)
  n=length(X.1)
  h0=p
  sampleVar=var(X.1)
  if (n<1) {
    stop("Invalid length of input data. It should be >0.0")
  }
  #X.1=as.data.frame(X.1)

  D4=sampleVar*((h0-0.5)^2*(3*(1-b-b^2)+b^3+b^4+b^5)+0.25*(1-(b+b^2+b^3+b^4+b^5)/5))
  SEW3=sqrt(abs((1/n)*(D4-sampleVar^2*((n-3)/(n-1)))))

  return(SEW3)

}
