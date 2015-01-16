#' Used to run the imageJ analysis using the tcltk interface

#' @description This function is used to run the imageJ analysis component of diskImageR using tcltk interface (i.e., pop-up boxes) to identify the project directory (the directory where all subsequent analyses will be saved) as well as the location of the directory containing the photos to be analyzed.

#' @param projectName the short name you want use for the project.
#' @param diskDiam the diameter of the diffusion disk in mm, defaults to 6.
#' @param imageJLoc the absolute path to imageJ on your computer. Current options are those standard for a mac: \code{imageJLoc} = "default" when imageJ is located at /Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub; \code{imageJLoc} = "loc2" for path /Applications/ImageJ.app/Contents/MacOS/JavaApplicationStub; \code{imageJ} will also accept any absolute path

#' @details Each photograph in the directory specified by \code{pictureDir} is input into ImageJ, where the built-in ‘find particles’ macro is used to find the center of a drug diffusion disk of the size specified by \code{diskDiam}. Lines are drawn every 5° out from the center of the disk, and the pixel intensity, which corresponds to cell density, is measured using the ‘plot-profile’ macro along each line. The results from all lines are saved into the "imageJ-out" directory in the specified \code{projectDir}. The average pixel intensity is then determined across all 72 lines at each distance and saved to \code{projectName}. \cr Note that the photograph names can be fairly important downstream and should follow a fairly strict convention to be able to take advantage of some of the built-in functions. Photographs should be named "line_factor1_factor2_factor3_...".

#' @return A .csv file is saved to the directory "imageJ-out" in the directory specified by \code{projectDir}. The average line for each photograph is saved to the list \code{projectName} in the global environment.

#' @export

#' @seealso \code{\link{runIJManual}} to run the imageJ analysis by manually supplying the project folder and photo location paths 

#'@author Aleeza C. Gerstein

runIJ <-
function(projectName, imageJLoc="default", diskDiam = 6){
	projectDir <- tk_choose.dir(caption = "Select main project directory") 
	inputDir <- tk_choose.dir(caption = "Select location of photographs")
#fill this in!
	# script <- .libPaths()/diskImageR/scr/
	script <- file.path(.libPaths(), "diskImageR/IJ_diskImageR.txt")
	# script <- tk_choose.files(caption = "Select the imageJ disk assay script")
	fileDir <- paste(Sys.Date(), projectName, sep="_")
	outputDir <- paste(projectDir, "/imageJ-out/", fileDir, "/", sep="")
	inputDir2 <- paste(inputDir, "/", sep="")
	IJarguments <- paste(inputDir2, outputDir, diskDiam, sep="*")

	if(length(dir(outputDir)) > 0){
		stop("Output files already exist in specified directory. Please delete existing files or change project name before continuing.")
	}
	
	dir.create(paste(projectDir, "/imageJ-out/", sep=""), showWarnings=FALSE)
	dir.create(paste(outputDir, sep=""), showWarnings= FALSE)
	dir.create(paste(projectDir, "/figures/", sep=""), showWarnings= FALSE)
	dir.create(paste(projectDir, "/parameter_files/", sep=""), showWarnings=FALSE)
			
	if (imageJLoc=="default" | imageJLoc=="loc2" ){
		if (imageJLoc=="loc2"){
			call <- paste("/Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")}
		if (imageJLoc=="default"){
			call <- paste("/Applications/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")}
	}
	else {call <- paste(imageJLoc,  "-batch", script, IJarguments, sep=" ")
		}
		
	system(call)

	cat(paste("\nOutput saved in directory: ", fileDir, "\n", sep=""))
	cat(paste("\nElements in dataframe ", projectName, ": \n", sep=""))	
	temp <- .ReadIn_DirCreate(projectDir, outputDir, projectName)
	cat("\a")
	assign(projectName, temp, envir=globalenv())
	}
	
# Create required directories in the main project directory and read in the imageJ output to a single list
# Helper function used to create new directories in the main project directory and read in the imageJ output created for each image to be stored in a single list. Used in \code{\link{runIJ}} and \code{\link{runIJManual}}

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
	newdir1 <- paste(workingDir, "/figures/", sep="")
	newdir2 <- paste(workingDir, "/parameter_files/", sep="")		
	if (!file.exists(newdir1)){
		dir.create(newdir1, showWarnings = FALSE)
		cat(paste("new directory: ", newdir1), sep="")
		}
	if (!file.exists(newdir2)){		
		dir.create(newdir2, showWarnings = FALSE)
		cat(paste("new directory: ", newdir3), sep="")
		}
	tList
	}

# Used to recursively read in and average imageJ output for each picture
#description Helper function used in \code{\link{.Readin_DirCreate}} to recursively read in the plot profile information for all 72 lines for each picture, find the average at each position, and save the average to a list 

.readIn <-
function(directoryPath, newList = list(), numDig=30) {
	currDir <- getwd()
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

# Used to load the imageJ output
# Helper function to load imageJ output in \code{\link{.readIn}}

.load.data <-
function(filename) {
	d <- read.csv(filename, header=TRUE, sep="\t")
   names(d) <- c("count", "distance","x")
   d
 }

