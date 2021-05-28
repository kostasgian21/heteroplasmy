#' A Function to read mouse heteroplasmy data (not finished!)
#'
#' This function allows you to read mouse heteroplasmy data from external files. Use with caution for now.
#' @param nameD Either "HB" or "LE".
#' @keywords heteroplasmy data
#' @export
#' @examples
#' readHeteroplasmyData(nameD="LE")

readHeteroplasmyData <- function(nameD="HB") {
  # return a new vector containing the mouse data
  if (nameD=="HB") {
    mouseData <- read.table("HB oocyte data.txt", sep="\t",header=T)
    mouseData=mouseData[,-c(1,2)]
    mouseData=mouseData[-c(1:11),]
    mouseData=mouseData[which(mouseData[,]!=""),]
    mouseData[, ] <- sapply(mouseData[, ], as.numeric)
  }else if(nameD=="LE"){
    mouseData <- read.table("LE oocyte data.txt", sep="\t",header=T)
    mouseData=mouseData[,-1]
    mouseData=mouseData[-c(1:10),]
    mouseData=mouseData[which(mouseData[,]!=""),]
    mouseData[, ] <- sapply(mouseData[, ], as.numeric)
  }else{
    mouseData=NULL
    print("wrong input, please use HB or LE")
  }
  return(mouseData)
}
