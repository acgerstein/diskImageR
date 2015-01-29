#' Run the imageJ analysis macro 

#' @description Used to run the imageJ analysis component of diskImageR and then load in the output from imageJ into R.

#' @param projectName the short name you want use for the project.
#' @param projectDir the path to the project directory where all analyses will be saved. If using tcltk interface to select directories leave as is (default=NA)
#' @param pictureDir the path to the directory where the pictures are to be analyzed. If using tcltk interface to select directories leave as is (default=NA)
#' @param diskDiam the diameter of the diffusion disk in mm, defaults to 6.
#' @param imageJLoc the absolute path to imageJ on your computer. Current options are those standard for a mac: \code{imageJLoc} = "default" when imageJ is located at /Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub; \code{imageJLoc} = "loc2" for path /Applications/ImageJ.app/Contents/MacOS/JavaApplicationStub; \code{imageJ} will also accept any absolute path. If you are on a windows machine please save imageJ to the default location when installing (Program Files/ImageJ).
#' @param manual whether to supply path names manually (defalt = TRUE) or using tcltk interface, i.e., pop-up boxes (manual = FALSE).


#' @details Each photograph in the directory specified by \code{pictureDir} is input into ImageJ, where the built-in 'find particles' macro is used to find the center of a drug diffusion disk of the size specified by \code{diskDiam}. Lines are drawn every 5 degrees out from the center of the disk, and the pixel intensity, which corresponds to cell density, is measured using the 'plot-profile' macro along each line. The results from all lines are saved into the "imageJ-out" directory in the specified \code{projectDir}. The average pixel intensity is then determined across all 72 lines at each distance and saved to \code{projectName}. \cr Note that the photograph names can be fairly important downstream and should follow a fairly strict convention to be able to take advantage of some of the built-in functions. Photographs should be named "line_factor1_factor2_factor3_...".

#' @return A .csv file is saved to the directory "imageJ-out" in the directory specified by \code{projectDir}. The average line for each photograph is saved to the list \code{projectName} in the global environment.

#' @export

#' @seealso \code{\link{runIJManual}} to run the imageJ analysis by manually supplying the project folder and photo location paths 

#'@author Aleeza C. Gerstein

IJMacro <-
function(projectName, projectDir=NA, pictureDir=NA, imageJLoc="loc2", diskDiam = 6, manual=TRUE){
	fileDir <- projectName
	if(is.na(projectDir)){
		projectDir <- tcltk::tk_choose.dir(caption = "Select main project directory") 
	}		
	if(is.na(pictureDir)){
		pictureDir <- tcltk::tk_choose.dir(caption = "Select location of photographs")
		pictureDir <- file.path(pictureDir, "")
	}
	
	dir.create(file.path(projectDir, "imageJ_out"), showWarnings=FALSE)
	outputDir <- file.path(projectDir, "imageJ_out", fileDir, "")
	IJarguments <- paste(pictureDir, outputDir, diskDiam, sep="*")	

	if(length(dir(outputDir)) > 0){
		cont <- readline(paste("Output files exist in directory ", outputDir, "\nOverwrite? (y/n) ", sep=""))
		if(cont=="n"){
			stop("Please delete existing files or change project name before continuing.")
			}
		if(cont=="y"){
			unlink(outputDir, recursive = TRUE)
		}
	}
	
	dir.create(file.path(outputDir), showWarnings= FALSE)
	dir.create(file.path(projectDir, "figures"), showWarnings=FALSE)
	dir.create(file.path(projectDir, "figures", fileDir), showWarnings=FALSE)
	dir.create(file.path(projectDir, "parameter_files"), showWarnings=FALSE)
	dir.create(file.path(projectDir, "parameter_files", fileDir), showWarnings=FALSE)	
			
	if(.Platform$OS.type=="windows"){
		script <- file.path(.libPaths(), "diskImageR", "IJ_diskImageR.ijm")[1]
		script <- gsub("Program Files", "progra~1", script)
		cmd <- "C:/progra~1/ImageJ/ij.jar"
		args <- paste("-batch", script, IJarguments)
		shell(paste(cmd, args))

	}
	else{		
		script <- file.path(.libPaths(), "diskImageR", "IJ_diskImageR.ijm")[1]
		if (imageJLoc=="default" | imageJLoc=="loc2" ){
			if (imageJLoc=="loc2"){
				call <- paste("/Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")}
			if (imageJLoc=="default"){
				call <- paste("/Applications/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")}
		}
		else {call <- paste(imageJLoc,  "-batch", script, IJarguments, sep=" ")
		}
		system(call)
	}
	
	cat(paste("\nOutput of imageJ analyses saved in directory: ", outputDir, "\n", sep=""))
	cat(paste("\nElements in dataframe ", projectName, ": \n", sep=""))	
	temp <- .ReadIn_DirCreate(projectDir, outputDir, projectName)
	cat("\a")
	assign(projectName, temp, envir=globalenv())
	}

.ReadIn_DirCreate <-
function(workingDir, folderLoc, experAbbr){
    setwd(workingDir)
	tList <- list()
	tList <- .readIn(folderLoc, tList, 30)
	len <- c()
		for (i in 1:length(tList)){
		len[i] <- length(tList[[i]][,1])
		}
	temp <- data.frame(names = names(tList), len)
	redo <- subset(temp, len==1, names)	
	tList
	}

.readIn <-function(directoryPath, newList = list(), numDig=30) {
	currDir <- getwd()
	print(currDir)
	getData <- function(i, newList, names) {
		if (i > length(dir())){
			names(newList) <- names
			print(names(newList))
			setwd(currDir)
			return (newList)
			}
		else {
			allLines <-  aggregate(.load.data(dir()[i])$x,  .load.data(dir()[i])["distance"], mean)
			newList[[length(newList)+1L]] <-  data.frame(distance = allLines[,1]*40/length(allLines[,1]), x= allLines[,2])
			temp <- paste(substr(basename(dir()[i]),1,numDig), "", sep="")
			names[i] <- strsplit(temp,".txt")[[1]][1]
			getData(i+1, newList, names)
		}
	}
	setwd(directoryPath)
	i <-1
	names <- c()
	findMin <- c()
	getData(i, newList, names)
}

.load.data <-
function(filename) {
	d <- read.csv(filename, header=TRUE, sep="\t")
   names(d) <- c("count", "distance","x")
   d
 }




