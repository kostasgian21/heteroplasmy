#' An example plotting function
#'
#' This function is used as a toy example on how to represent the data statistics
#' regarding the variance of the sample. The mean variance and its standard error
#' are depicted. Note that this is just an illustration to show that the analytic
#' and the resampling approaches almost match each other.
#' @param functions Choose the subset of the functions you wish use for the calculation and
#' subsequent plot of the standard error of the variance. You can use one or a combination
#' of "normalApr", "analytic","bootstrap","correctedBoot", and "jackknife". For now,
#' it outputs all of the aforementioned methods!
#' @param data The input data in the form of a dataframe or matrix (which will be transformed into
#' a dataframe).
#' @keywords plot,standard,error
#' @export
#' @examples
#' # size of the sample
#' n=50
#' #generate a random sample of size n from a normal distribution
#' data_ex=rnorm(n,0.5,0.1)
#' plotStdErrVar(data_ex)

plotStdErrVar <- function(data,functions=c("normalApr","analytic","bootstrap","correctedBoot","jackknife"),...) {
  functions <- list(...)
  if (length(functions)==0) {
    functions=c("normalApr","analytic","bootstrap","correctedBoot","jackknife")
  }
  sampleVar=var(data)

  resultMeans= data.frame(W1=numeric(),
                          W2=numeric(),
                          Boots=numeric(),
                          CorBoots=numeric(),
                          Jack=numeric(),
                          SEW1=numeric(),
                          SEW2=numeric(),
                          SEBoots=numeric(),
                          SECorBoots=numeric(),
                          SEJack=numeric())

  if ("analytic" %in% functions) {
    anvar=analyticVar(data)
  }
  if ("normalApr" %in% functions) {
    anvarN=analyticVar(data,normal = TRUE)
  }
  if ("bootstrap" %in% functions) {
    bootVar=bootstrapVar(data)
  }
  if ("correctedBoot" %in% functions) {
    bootVarOpt=bootstrapVar(data,sigmaOpt = TRUE)
  }

  if ("jackknife" %in% functions) {
    jackVar=jackVar(data)
  }
  resultMeans[1,] = c(sampleVar,sampleVar,sampleVar,sampleVar,
                      sampleVar,anvarN,anvar,bootVar,
                      bootVarOpt,jackVar)

  ns=c(1)
  # row.names(resultMeans) <- ns
  cl <- c("red","forestgreen","blue","gold","brown")#rainbow(4)
  # cl <- rainbow(5)
  names(resultMeans)[1] <- "W1 normal appr"
  names(resultMeans)[2] <- "Wonnapinij et al."
  names(resultMeans)[3] <- "corBoots"
  names(resultMeans)[4]<-"corBoots with g(x)"

  #plot the results
  for (i in 1:5){
    plot(y=resultMeans[,i],x=ns,ylab="Var(h)",xlab="#instance",
         col = cl[i],type = "p",main="h from HB oocyte data",cex =1.2,
         # col = cl[i],type = "p",main=TeX(r'(h from HB oocyte data $\mu$)'),cex =1.2,
         ylim = c(0,max(resultMeans[,1:5])*2))
    for (k in 1:length(ns)) {
      arrows(x0=ns[k],x1=ns[k], y0=resultMeans[k,i]-resultMeans[k,i+5], y1=resultMeans[k,i]+resultMeans[k,i+5],
             col = cl[i],code=3, angle=90, length=0.1)
    }
    par(new=TRUE)
  }
  legend("topright", legend=colnames(resultMeans[,1:5]),
         col=cl, lty=1:2, cex=0.8)#box.lty=0
}
