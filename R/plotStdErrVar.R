#' A Function to read mouse heteroplasmy data
#'
#' This function allows you to read mouse heteroplasmy data from external files.
#' @param nameD Either "HB" or "LE".
#' @keywords heteroplasmy data
#' @export
#' @examples
#' plotStdErrVar(data)

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
  for (i in 1:4){
    plot(y=resultMeans[,i],x=ns,ylab="Var(h)",xlab="#instance",
         col = cl[i],type = "p",main="h from HB oocyte data",cex =1.2,
         # col = cl[i],type = "p",main=TeX(r'(h from HB oocyte data $\mu$)'),cex =1.2,
         ylim = c(0,max(resultMeans[,1:4])*2))
    for (k in 1:length(ns)) {
      arrows(x0=ns[k],x1=ns[k], y0=resultMeans[k,i]-resultMeans[k,i+5], y1=resultMeans[k,i]+resultMeans[k,i+5],
             col = cl[i],code=3, angle=90, length=0.1)
    }
    par(new=TRUE)
  }
  legend("topright", legend=colnames(resultMeans[,1:4]),
         col=cl, lty=1:2, cex=0.8)#box.lty=0

}
