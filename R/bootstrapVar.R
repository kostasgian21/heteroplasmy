#' A Function to read mouse heteroplasmy data
#'
#' This function allows you to read mouse heteroplasmy data from external files.
#' @param nameD Either "HB" or "LE".
#' @keywords heteroplasmy data
#' @export
#' @examples
#' bootstrapVar(data)

bootstrapVar <- function(data,nrep=1000,biased=FALSE,corrected=FALSE,sigmaOpt=FALSE) {
  n=length(which(!is.na(data)))
  X.1 = data
  #X.1=X.1/100
  h0=mean(X.1)
  opti = 2.61
  coeB=-4.1434
  coeC=-0.7387
  coeD=1.0125
  coeE=2.8390

  X.1_boot=as.data.frame(X.1)
  boots = vector()
  bootVars=vector()
  bootVarsBiased=vector()

  for (i in 1:nrep) {
    repeat {
      boot= sample(1:nrow(X.1_boot), nrow(X.1_boot), replace=TRUE)
      boots=X.1_boot[boot,1]
      if (!is.na(var(boots,na.rm = TRUE))) {
        break
      }
    }
    #append the corrected variance of the i-th bootstrap sample
    bootVars= c(bootVars,var(boots,na.rm = TRUE)*(n/(n - 1)))
    bootVarsBiased= c(bootVarsBiased,var(boots,na.rm = TRUE))

  }

  # standard error of the bootstrap samples
  SEB=sd(bootVars)
  heteropVar=SEB
  if (biased==TRUE) {
    warning("Cautious, biased calculation always underestimates the variance!")
    SEB=sd(bootVarsBiased)
    heteropVar=SEB
  }

  if (corrected==TRUE) {
    message("Correceted by a constant factor")
    heteropVar=SEB*((n - opti)/n)
  }
  if (sigmaOpt==TRUE) {
    message("Correceted by a fitted sigmoid function, ie a four-parameter log-logistic function")
    heteropVar=SEB*((n - opti)/n)/(coeC+(coeD-coeC)/(1+exp(coeB*(log(n)-log(coeE)))))

  }

  return(heteropVar)
}
