#' A Function to read mouse heteroplasmy data (not finished!)
#'
#' This function allows you to read mouse heteroplasmy data from external files.
#' Use with caution (for now). Calling this function is equivalent to loading the already stored dataframe with heteroplasmy measurements (ie, \emph{mousedataHB},\emph{mousedataLE},\emph{mousedataFreyer}).
#' @param nameD Either "HB" or "LE" or "Freyer" or "Broz".
#' @return A dataframe containing mouse heteroplasmy data.
#' @keywords heteroplasmy data
#' @export
#' @examples
#' df1=readHeteroplasmyData(nameD="LE")
#' df2=readHeteroplasmyData(nameD="HB")
#' df3=readHeteroplasmyData(nameD="Freyer")
#' df4=readHeteroplasmyData(nameD="Broz")
#'
#' # or equivalent
#' df1=mousedataLE
#' df2=mousedataHB
#' df3=mousedataFreyer
#' df4=mousedataBroz

readHeteroplasmyData <- function(nameD="HB") {
  # return a new vector containing the mouse data
  if (nameD=="HB") {
    mouseData <- read.table("HB oocyte data.txt", sep="\t",header=T)
    mouseData=mouseData[,-c(1,2)]
    mouseData=mouseData[-c(1:11),]
    mouseData=mouseData[which(mouseData[,]!=""),]
    mouseData[, ] <- sapply(mouseData[, ], as.numeric)
    mouseData=mouseData[1:25,]
    rnames=seq(1,25)
    row.names(mouseData) <- rnames
    mouseData=mouseData/100
  }else if(nameD=="LE"){
    mouseData <- read.table("LE oocyte data.txt", sep="\t",header=T)
    mouseData=mouseData[,-1]
    mouseData=mouseData[-c(1:10),]
    mouseData=mouseData[which(mouseData[,]!=""),]
    mouseData[, ] <- sapply(mouseData[, ], as.numeric)
    mouseData=mouseData[1:20,]
    rnames=seq(1,20)
    row.names(mouseData) <- rnames
    mouseData=mouseData/100
  }else if(nameD=="Freyer"){
    mouseData <- read.table("freyerPGC.txt", sep = "\t",
                            header = F)
  }else if(nameD=="Broz"){
    mouseData <- read.table("orgDatBroz.txt", sep = "\t",
                            header = T)
  }else{
    mouseData=NULL
    stop("wrong input, please use HB or LE")
  }
  return(mouseData)
}
