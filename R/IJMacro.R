#' Run the imageJ analysis macro 

#' @description Used to run the imageJ analysis component of diskImageR and then load in the output from imageJ into R.

#' @param projectName the short name you want use for the project.
#' @param projectDir the path to the project directory where all analyses will be saved. If using tcltk interface to select directories leave as is (default=NA)
#' @param photoDir the path to the directory where the photographs are to be analyzed. If using tcltk interface to select directories leave as is (default=NA)
#' @param diskDiam the diameter of the diffusion disk in mm, defaults to 6.
#' @param imageJLoc the absolute path to imageJ on your computer. Current options are those standard for a mac: \code{imageJLoc} = "default" when imageJ is located at /Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub; \code{imageJLoc} = "loc2" for path /Applications/ImageJ.app/Contents/MacOS/JavaApplicationStub. If you are on a windows machine please save imageJ to the default location when installing (Program Files/ImageJ). If you wish to run imageJ from an alternative path use \code{imageJLoc} to specify the absolute path.

#' @details Each photograph in the directory specified by \code{photoDir} is input into ImageJ, where the built-in 'find particles' macro is used to find the center of a drug diffusion disk of the size specified by \code{diskDiam}. Lines are drawn every 5 degrees out from the center of the disk, and the pixel intensity, which corresponds to cell density, is measured using the 'plot-profile' macro along each line. The results from all lines are saved into the "imageJ-out" directory in the specified \code{projectDir}. The average pixel intensity is then determined across all 72 lines at each distance and saved to \code{projectName}. \cr Note that the photograph names can be fairly important downstream and should follow a fairly strict convention to be able to take advantage of some of the built-in functions. Photographs should be named "line_factor1_factor2_factor3_...".

#' @return A .csv file is saved to the directory "imageJ-out" in the directory specified by \code{projectDir}. The average line for each photograph is saved to the list \code{projectName} in the global environment.

#' @export

#'@author Aleeza C. Gerstein

IJMacro <-
function(projectName, projectDir=NA, photoDir=NA, imageJLoc="loc2", diskDiam = 6){
	# if(!is.char(projectName))
	fileDir <- projectName
	if(is.na(projectDir)){
		projectDir <- tcltk::tk_choose.dir(caption = "Select main project directory") 
	}		
	if(is.na(photoDir)){
		photoDir <- tcltk::tk_choose.dir(caption = "Select location of photographs")
		photoDirOrig <- photoDir
		photoDir <- file.path(photoDir, "")
		if (projectDir == photoDirOrig) {
			cat("The photograph directory can not be used for the main project directory. Please select a different folder for the main project directory.")
			projectDir <- tcltk::tk_choose.dir(caption = "Select main project directory") 
		}

	}
	setwd(photoDir)
	if (TRUE %in% file.info(dir())[,2]) {
		stop("There is a folder located in your photograph directory. Please remove before continuing.")
		}
	dir.create(file.path(projectDir, "imageJ_out"), showWarnings=FALSE)
	outputDir <- file.path(projectDir, "imageJ_out", fileDir, "")
	IJarguments <- paste(photoDir, outputDir, diskDiam, sep="*")	

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
	
	script <- file.path(.libPaths(), "diskImageR", "IJ_diskImageR.ijm")[1]			
	if(.Platform$OS.type=="windows"){
		IJarguments <- paste(paste(photoDir,  "", sep="\\"), paste(outputDir, "", sep="\\"), diskDiam, sep="*")		
		script <- gsub("Program Files", "progra~1", script)
		knownIJLoc <- FALSE
		if("ImageJ.exe" %in% "C:\\progra~1\\ImageJ\\"){
		  	cmd <- "C:\\progra~1\\ImageJ\\ImageJ.exe"
		  	knownIJLoc <- TRUE
		  	}
		if("ImageJ.exe" %in% "C:\\Program Files (x86)\\ImageJ\\"){
			cmd <- '"C:\\Program Files (x86)\\ImageJ\\ImageJ.exe"'
			knownIJLoc <- TRUE
			}
		if("ImageJ.exe" %in% imageJLoc){
			cmd <- paste(imageJLoc, "ImageJ.exe", sep="")
			knownIJLoc <- TRUE
			}
		if(knownIJloc == FALSE){
			stop("ImageJ is not in expected location. Please move ImageJ to the Program Files directory, or specify the path to its location using the argument 'imageJLoc'")
		}		
		args <- paste("-batch", script, IJarguments)
		args <- gsub("/", "\\\\", args)
		shell(paste(cmd, args), wait=TRUE,intern=TRUE)
	}
	else{				
		if (imageJLoc=="default" | imageJLoc=="loc2" ){
			if ("ImageJ.app" %in% dir("/Applications/")){
				call <- paste("/Applications/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")}
			if ("ImageJ.app" %in% dir("/Applications/ImageJ/")){			
				call <- paste("/Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")}
		}
		else {
			if ("ImageJ.app" %in% imageJLoc){
				call <- paste(imageJLoc,  "-batch", script, IJarguments, sep=" ")
				}
			else{
			stop("ImageJ is not in expected location. Please move ImageJ to the Applications directory, or specify the path to its location using the argument 'imageJLoc'")
				}
	}
		system(call)
		}

	count_wait<-0.0;
	while(length(dir(outputDir))<length(dir(photoDir)) && count_wait<1e12)
	{
	  count_wait<-count_wait+1.0
	}
	cat(paste("\nOutput of imageJ analyses saved in directory: \n", outputDir, "\n", sep=""))
	cat(paste("\nElements in list '", projectName, "': \n", sep=""))	
	temp <- .ReadIn_DirCreate(projectDir, outputDir, projectName)
	if(!length(dir(photoDir)) == length(temp)){
		stop("Mismatch between the number of files in the photograph directory and the number of images analyzed. This likely indicates a non-photograph file is located in this directory. Please remove and rerun before continuing.")
		}		
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
	# print(currDir)
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




