#' Used to run the imageJ analysis by manually supplying path locations

#' @description This function is used to run the imageJ analysis component of diskImageR by manually supplying the path location of the project directory (the directory where all subsequent analyses will be saved) and the location of the directory containing the photos to be analyzed.

#' @param projectName the short name you want use for the project.
#' @param projectDir the path to the project directory where all analyses will be saved
#' @param pictureDir the path to the directory where the pictures are to be analyzed
#' @param diskDiam the diameter of the diffusion disk in mm, defaults to 6.
#' @param imageJLoc the absolute path to imageJ on your computer. Current options are standard for a mac: \code{imageJLoc} = "default" when imageJ is located at /Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub; ly \code{imageJLoc} = "loc2" for path /Applications/ImageJ.app/Contents/MacOS/JavaApplicationStub

#' @details Each photograph in the directory specified by \code{pictureDir} is input into ImageJ, where the built-in 'find particles' macro is used to find the center of a drug diffusion disk of the size specified by \code{diskDiam}. Lines are drawn every 5 degrees out from the center of the disk, and the pixel intensity, which corresponds to cell density, is measured using the 'plot-profile' macro along each line. The results from all lines are saved into the "imageJ-out" directory in the specified \code{projectDir}. The average pixel intensity is then determined across all 72 lines at each distance and saved to \code{projectName}.

#' @return A .csv file is saved to the directory "imageJ-out" in the directory specified by \code{projectDir}. The average line for each photograph is saved to the list \code{projectName} in the global environment.

#' @export

#' @seealso \code{\link{runIJ}} to run the imageJ analysis using the tcltk interface (i.e., pop-up boxes) to supply the project and picture directories. Test

#'@author Aleeza C. Gerstein

runIJManual <-
function(projectName, projectDir, pictureDir, imageJLoc="loc2", diskDiam = 6){
	fileDir <- projectName
	outputDir <- file.path(projectDir, "imageJ-out", fileDir, "")
	script <- file.path(.libPaths(), "diskImageR", "IJ_diskImageR.txt")
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
	
	
	dir.create(file.path(projectDir, "imageJ-out"), showWarnings=FALSE)
	dir.create(outputDir, showWarnings= FALSE)
	dir.create(file.path(projectDir, "figure"), showWarnings= FALSE)
		
	if (imageJLoc=="default" | imageJLoc=="loc2" ){
		if (imageJLoc=="loc2"){
			call <- paste("/Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")}
		if (imageJLoc=="default"){
			call <- paste("/Applications/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")}
	}
	else {call <- paste(imageJLoc,  "-batch", script, IJarguments, sep=" ")
		}
	print(call)
	system(call)

	cat(paste("\nOutput of imageJ analyses saved in directory: ", outputDir, "\n", sep=""))
	cat(paste("\nElements in dataframe ", projectName, ": \n", sep=""))	
	temp <- .ReadIn_DirCreate(projectDir, outputDir, projectName)
	cat("\a")
	assign(projectName, temp, envir=globalenv())
	}
	
