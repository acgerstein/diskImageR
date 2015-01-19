#' Read in an existing dataframe using the tcltk interface

#' @description Open an existing dataframe previously created with either \code{createDataframe} or \code{aggregateData} using tcltk interface.

#' @param projectName the short name you want use for the project. Note, this can be different than the previous file name used for this project.

#' @return \code{projectName} is saved to the global directory and can be used for \code{aggregateData}

#' @export

#' @seealso \code{\link{runIJManual}} to run the imageJ analysis by manually supplying the project folder and photo location paths 

#'@author Aleeza C. Gerstein

readExistingDF <- function(projectName){
	projectFolder <- tcltk::tk_choose.dir(caption = "Select main project directory")
	setwd(projectFolder)
	fileName <- tcltk::tclvalue(tcltk::tkgetOpenFile())
	if (!nchar(fileName)){
		tcltk::tkmessageBox(message = "No file was selected")
	}
	newDF <- read.csv(fileName)
	if (grepl(".df", fileName) == TRUE){
		assign(paste(projectName, ".df", sep=""), newDF, envir=globalenv())
		name <- paste(projectName, ".df", sep="")
	}
	else {
		assign(paste(projectName, ".ag", sep=""), newDF, envir=globalenv())	
		name <- paste(projectName, ".ag", sep="")
		}

	cat(paste(name, "now present in working directory"))
}