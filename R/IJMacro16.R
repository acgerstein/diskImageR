#' Run an imageJ analysis macro on the folder that contains the photograph to be analyzed

#' @description \code{IJMacro} is used to run the imageJ analysis component of diskImageR and then load in the acquired output from imageJ into R.

#' @param projectName the short name you want use for the project
#' @param projectDir the path to the project directory where all analyses will be saved. If left as NA (the default) you will be able to specify the locaion through a pop-up box. (default=NA)
#' @param photoDir the path to the directory where the photographs are to be analyzed. If left as NA (the default) you will be able to specify the locaion through a pop-up box. (default=NA)
#' @param diskDiam the diameter of the diffusion disk in mm, defaults to 6.
#' @param imageJLoc the absolute path to ImageJ (\href{http://rsb.info.nih.gov/ij/download.html}{ImageJ}) on your computer. Leave as NA (the default) if you have downloaded ImageJ to a standard location (Mac: /Applications/ImageJ.app or /Applications/ImageJ/ImageJ.app/; Windows: Program Files/ImageJ). If you wish to run imageJ from an alternative path use \code{imageJLoc} to specify the absolute path.

#' @details Each photograph in the directory specified by \code{photoDir} is input into ImageJ, where the built-in 'find particles' macro is used to find the center of a drug diffusion disk of the size specified by \code{diskDiam}. Lines are drawn every 5 degrees out from the center of the disk, and the pixel intensity, which corresponds to cell density, is measured using the 'plot-profile' macro along each line. The results from all lines are saved into the "imageJ-out" directory in the specified \code{projectDir}. The average pixel intensity is then determined across all 72 lines for each photograph and saved to \code{projectName}. \cr Note that the photograph names can be fairly important downstream and should follow a fairly strict convention to be able to take advantage of some of the built-in functions. Photographs should be named "line_factor1_factor2_factor3_...".

#' @section Important:
#' There can not be any spaces or special characters in any of the folder names that are in the path that lead to either the main project directory or the photograph directory. If there are an error box titled "Macro Error" will pop up and the script will not run.
#' The project name should ideally be fairly short (easy to type without typos!) and specific to the project. It must start with a letter, not a number or special character, but can otherwise be anything. The project name must always be specified with quotation marks around it (a surprisingly common error).

#' @return A .csv file is saved to the directory "imageJ_out" in the directory specified by \code{projectDir}. The average line for each photograph is saved to the list \code{projectName} in the global environment.

#' @examples
#' \dontrun{
#' IJMacro("myProject")
#' }

#' @export

IJMacro16 <-
function(projectName, projectDir=NA, photoDir=NA, imageJLoc=NA, diskDiam = 6, drugs = c("CIP5", "S10", "FOX30", "NA30", "CT10", "IPM10", "AMC30", "CTX30", "F300", "TE30", "CAZ30", "C30", "SXT25", "CPD10", "AM10", "ATM30")){
	# if(!is.char(projectName))
	diskImageREnv <- new.env()
	fileDir <- projectName
	#get R version and use appropriate tcltk
	if(is.na(projectDir)){
		projectDir <- tcltk::tk_choose.dir(caption = "Select main project directory")
		if(is.na(projectDir)) stop("")
		}
	if(is.na(photoDir)){
		photoDir <- tcltk::tk_choose.dir(caption = "Select location of photographs")
		if(is.na(photoDir)) stop("")
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
	dir.create(file.path(projectDir, "disk_coordinates"), showWarnings=FALSE)
	outputDir <- file.path(projectDir, "imageJ_out", fileDir, "")
	outputDirCoord <- file.path(projectDir, "disk_coordinates", fileDir, "")
	IJarguments <- paste(photoDir, outputDir, outputDirCoord, diskDiam, sep="*")

	if(length(dir(outputDir)) > 0){
		cont <- readline(paste("Output files exist in directory ", outputDir, "\nOverwrite? [y/n] ", sep=""))
		if(cont=="n"){
			stop("Please delete existing files or change project name before continuing.")
			}
		if(cont=="y"){
			unlink(outputDir, recursive = TRUE)
		}
	}

	dir.create(file.path(outputDir), showWarnings= FALSE)
	dir.create(file.path(outputDirCoord), showWarnings= FALSE)
	dir.create(file.path(projectDir, "figures"), showWarnings=FALSE)
	dir.create(file.path(projectDir, "figures", fileDir), showWarnings=FALSE)
	dir.create(file.path(projectDir, "parameter_files"), showWarnings=FALSE)
	dir.create(file.path(projectDir, "parameter_files", fileDir), showWarnings=FALSE)

	script <- file.path(.libPaths(), "diskImageR", "IJ_diskImageR16.ijm")[1]
	if(.Platform$OS.type=="windows"){
		IJarguments <- paste(paste(photoDir,  "", sep="\\"), paste(outputDir, "", sep="\\"),paste(outputDirCoord, "", sep="\\"),  diskDiam, sep="*")
		script <- gsub("Program Files", "progra~1", script)
		knownIJLoc <- FALSE
		if("ImageJ.exe" %in% dir("C:\\progra~1\\ImageJ\\")){
		  	cmd <- "C:\\progra~1\\ImageJ\\ImageJ.exe"
		  	knownIJLoc <- TRUE
		  	}
		if("ImageJ.exe" %in% dir("C:\\Program Files (x86)\\ImageJ\\")){
			cmd <- '"C:\\Program Files (x86)\\ImageJ\\ImageJ.exe"'
			knownIJLoc <- TRUE
			}
		if("ImageJ.exe" %in% imageJLoc){
			cmd <- paste(imageJLoc, "ImageJ.exe", sep="")
			knownIJLoc <- TRUE
			}
		if(knownIJLoc == FALSE){
			stop("ImageJ is not in expected location. Please move ImageJ to the Program Files directory, or specify the path to its location using the argument 'imageJLoc'")
		}
		args <- paste("-batch", script, IJarguments)
		args <- gsub("/", "\\\\", args)
		shell(paste(cmd, args), wait=TRUE,intern=TRUE)
	}
	else{
		knownIJLoc <- FALSE
		if ("ImageJ.app" %in% dir("/Applications/")){
			call <- paste("/Applications/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")
			knownIJLoc <- TRUE
			}

		if (knownIJLoc == FALSE & "ImageJ.app" %in% dir("/Applications/ImageJ/")){
			call <- paste("/Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub -batch", script, IJarguments, sep=" ")
			knownIJLoc <- TRUE
			}
		if (knownIJLoc == FALSE & "ImageJ.app" %in% imageJLoc){
				call <- paste(imageJLoc,  "-batch", script, IJarguments, sep=" ")
				knownIJLoc <- TRUE
				}
		if(knownIJLoc == FALSE){
			stop("ImageJ is not in expected location. Please move ImageJ to the Applications directory, or specify the path to its location using the argument 'imageJLoc'")
				}
		system(call)
		}

	count_wait<-0.0;
	while(length(dir(outputDir))<length(dir(photoDir)) && count_wait<1e12)
	{
	  count_wait<-count_wait+1.0
	}

	cat("\n")
	cat("\n Processing data to determine which lines to use: \n")
	temp <- .ReadIn_DirCreate_Top(projectDir, outputDir, projectName)
	if(!length(dir(photoDir))*16 == length(temp)){
		stop("Mismatch between the number of files in the photograph directory and the number of images analyzed. This likely indicates a non-photograph file is located in this directory. Please remove and rerun before continuing.")
		}
	cat("\a")
	assign(projectName, temp, inherits=TRUE)
# cat("Assigning drug to disk coordinates")
.getCoordinates(projectName, drugs=drugs)
	 cat(paste("\nThe average line from each phogograph has been saved to: \n", file.path(getwd(), "parameter_files", projectName, paste("averageLines.csv", sep="")), "\n", sep=""))
	}

.saveAveLine <- function(L){
  addNA <- function(x, maxLen){
  	if(nrow(x) < maxLen){
  		diffLen <- maxLen - nrow(x)
	     tdf <- data.frame(rep(NA, diffLen), rep(NA, diffLen))
	     names(tdf) <- names(x)
  		x <- rbind(round(x, 3), tdf)
  	 }
  	else x <- round(x, 3)
  }
  maxLen <- max(sapply(L, nrow))
  newList <- lapply(L, addNA, maxLen)
   df <- data.frame(matrix(unlist(newList), nrow=maxLen))
   names(df) <- paste(c("distance", "instensity"), rep(names(L), each=2), sep=".")
   df
}

.ReadIn_DirCreate_Top <-
function(workingDir, folderLoc, experAbbr){
    setwd(workingDir)
	tList <- list()
	tList <- .readInTop(folderLoc, tList, 30)
	len <- c()
		for (i in 1:length(tList)){
		len[i] <- length(tList[[i]][,1])
		}
	temp <- data.frame(names = names(tList), len)
	redo <- subset(temp, len==1, names)
	tList
	}

.readInTop <-function(directoryPath, newList = list(), numDig=30, numTop = 20) {
	currDir <- getwd()
	# print(currDir)
	getData <- function(i, newList, names) {
		# print(i)
		if (i > length(dir())){
			names(newList) <- names
			setwd(currDir)
			return (newList)
			}
		else {
			cat(".")
			lines <-  data.frame(.load.data(dir()[i])$x,  .load.data(dir()[i])["distance"])
			names(lines) <- c("x", "distance")
			numPts <- length(lines$distance)/180
			newd <- data.frame(x = unique(lines$distance), L1 = lines$x[1:numPts])
		  start <- seq(1, length(lines$distance), by=length(lines$distance)/180)
		  for(j in 2:179){
		 	 newd <- cbind(newd, lines$x[start[j]:(start[j+1]-1)])
		  }
		 names(newd)[3:180] <- paste0("L",2:179)
		  # aveSorted <- apply(newd, 1, function(x) mean(sort(x)[(180-numTop):180]))
			# aveSorted <- apply(newd, 1, function(x) median(x))

			sum1 <- apply(newd[2:180], 2, function(x) sum(x[1:length(x)/2.5]))
			sum2 <- apply(newd[2:180], 2, function(x) sum(x[(length(x)/2.5): length(x)]))


			close <- which(sum1 > quantile(sum1, 0.4) & sum1 < quantile(sum1, 0.6))
			far <- which(sum2 > quantile(sum2, 0.75) & sum2 < quantile(sum2, 0.95))
			aveSorted_close <- apply(newd[1:(length(newd[,1])/2.5), close], 1, function(x) mean(x))
			aveSorted_far <- apply(newd[(length(newd[,1])/2.5):length(newd[,1]), far], 1, function(x) mean(x))
			# aveSorted <-  data.frame(distance = newd$x*28/length(newd$x), x= c(aveSorted_close, aveSorted_far))

			# ls0.8_close <- which(sum1 > quantile(sum1, 0.75) & sum1 < quantile(sum1, 0.85))
			# ls0.8_far <- which(sum2 > quantile(sum2, 0.75) & sum2 < quantile(sum2, 0.95))
			# aveSorted0.8_close <- apply(newd[1:(length(newd[,1])/2), ls0.8_close], 1, function(x) mean(x))
			# aveSorted0.8_far <- apply(newd[(length(newd[,1])/2+1):length(newd[,1]), ls0.8_far], 1, function(x) mean(x))
			# temp0.8b<-  data.frame(distance = newd$x*28/length(newd$x), x= c(aveSorted0.8_close, aveSorted0.8_far))

			# ls0.8 <- which(sum > quantile(sum, 0.7) & sum < quantile(sum, 0.9))
			# aveSorted0.8 <- apply(newd[, ls0.8], 1, function(x) mean(x))
			# temp0.8 <-  data.frame(distance = newd$x*28/length(newd$x), x= aveSorted0.8)
			# plot(temp0.8, ylim=c(0, 250), type="l", col="red")
			# points(temp0.8b, type="l")
      #
			# ls0.5_close <- which(sum1 > quantile(sum1, 0.4) & sum1 < quantile(sum1, 0.6))
			# ls0.5_far <- which(sum2 > quantile(sum2, 0.4) & sum2 < quantile(sum2, 0.6))
			# aveSorted0.5_close <- apply(newd[1:(length(newd[,1])/2), ls0.5_close], 1, function(x) mean(x))
			# aveSorted0.5_far <- apply(newd[(length(newd[,1])/2+1):length(newd[,1]), ls0.5_far], 1, function(x) mean(x))
			# temp0.5b<-  data.frame(distance = newd$x*28/length(newd$x), x= c(aveSorted0.5_close, aveSorted0.5_far))
      #
			# ls0.5 <- which(sum > quantile(sum, 0.4) & sum < quantile(sum, 0.6)) #35 lines
			# aveSorted0.5 <- apply(newd[, ls0.5], 1, function(x) mean(x))
			# temp0.5 <-  data.frame(distance = newd$x*28/length(newd$x), x= aveSorted0.5)
			# points(temp0.5, type="l", col="red")
			# points(temp0.5b, type="l", col="red")
      #
			# ls0.3 <- which(sum > quantile(sum, 0.2) & sum < quantile(sum, 0.3))
			# aveSorted0.3 <- apply(newd[, ls0.3], 1, function(x) mean(x))
			# temp0.3 <-  data.frame(distance = newd$x*28/length(newd$x), x= aveSorted0.3)
			# points(temp0.3)



			#the 28 comes from the IJ16 macro
			newList[[length(newList)+1L]] <-  data.frame(distance = newd$x*28/length(newd$x), x= c(aveSorted_close, aveSorted_far))
			# newList[[length(newList)+1L]] <-  data.frame(distance = newd$x*28/length(newd$x), x= aveSorted)
			# temp <- paste(substr(basename(dir()[i]),1,numDig), "", sep="")
			temp <- dir(directoryPath)[i]
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

# #testing what lines to keep (take out eventually)
# lines <- read.table("imageJ_out/disk16/ATCC25922_2.txt", sep="\t", header=TRUE, row.names=1)
# names(lines) <- c("distance", "x")
# numPts <- length(lines$distance)/180
# newd <- data.frame(x = unique(lines$distance), L1 = lines$x[1:numPts])
# start <- seq(1, length(lines$distance), by=length(lines$distance)/180)
# for(j in 2:179){
#  newd <- cbind(newd, lines$x[start[j]:(start[j+1]-1)])
# }
# names(newd)[3:180] <- paste0("L",2:179)
#
# par(mfrow=c(4, 2), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
# for(numTop in c(10, 20, 30, 40, 50, 60, 70, 80)){
# 	aveSorted <- apply(newd, 1, function(x) mean(sort(x)[(180-numTop):180]))
# 	plot(newd[,1], aveSorted, xaxt="n", yaxt="n", ylim=c(0, 250), main = paste0("numTop = ", numTop))
# 	abline(v=80)
# 	# plot(newd[,1]/7.9667, aveSorted, xaxt="n", yaxt="n", ylim=c(0, 250), main = paste0("numTop = ", numTop))
# 	if(numTop %in% c(10, 30, 50, 70)) axis(2, las=2)
# 	else axis(2, labels=FALSE)
# 	if(numTop > 60) axis(1)
# 	else axis(1, labels=FALSE)
# }
#
# sort(newd[80, 2:180])
# hist(newd[80, 2:180])

.load.data <-
function(filename) {
	d <- read.csv(filename, header=TRUE, sep="\t")
   names(d) <- c("count", "distance","x")
   d
 }

 .getCoordinates <- function(projectName, drugs){
   data <- eval(parse(text=projectName))
   d <- data.frame()
   d$line <- unlist(lapply(as.character(d$name), function(x) strsplit(x, "_")[[1]][2]))
   photoNames <- unique(unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1])))
    fileFolder <- projectName
    mapDir <- file.path(getwd(), "disk_coordinates", fileFolder)
		mapList <- list()
		photoNames <- unique(unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1])))
		i <- 0
		for (m in photoNames){
			i <- i +1
			mapList[[i]] <- read.csv(file.path(mapDir, paste0(m, "_ResultsTable.txt")), sep="\t")
			mapList[[i]]$photoName <- rep(m, 16)
			mapList[[i]]$XYpos <- c(order(mapList[[i]][1:4, "X"]), order(mapList[[i]][5:8, "X"])+4, order(mapList[[i]][9:12, "X" ])+8, order(mapList[[i]][13:16, "X" ])+12)
			mapList[[i]]$drug <- drugs
			mapList[[i]] <- mapList[[i]][order(mapList[[i]]$XYpos),]
		}
		 map <- do.call(rbind.data.frame, mapList)
     write.table(map, file.path(mapDir, paste0(projectName, "_ResultsTable.txt")), row.names=FALSE, sep="\t")
		 assign(paste(projectName, "map", sep="."), map, inherits=TRUE)
		 cat("\n")
		 cat(paste0("\n Drug coordinates map saved to: ", mapDir, projectName, "_ResultsTable.txt"))
		 cat(paste0("\n", projectName, ".map has been written to the global environment\n"))
     }


#  lines <-  data.frame(.load.data(dir()[i])$x,  .load.data(dir()[i])["distance"])
#  names(lines) <- c("x", "distance")
#  numPts <- length(lines$distance)/72
#  newd <- data.frame(x = unique(lines$distance), L1 = lines$x[1:numPts])
#  start <- seq(1, length(lines$distance), by=length(lines$distance)/72)
#  for(j in 2:71){
# 	newd <- cbind(newd, lines$x[start[j]:(start[j+1]-1)])
#  }
# names(newd)[3:72] <- paste0("L",2:71)
#
#  aveSorted <- apply(newd, 1, function(x) mean(sort(x)[62:72]))

# .readIn <-function(directoryPath, newList = list(), numDig=30) {
# 	currDir <- getwd()
# 	# print(currDir)
# 	getData <- function(i, newList, names) {
# 		if (i > length(dir())){
# 			names(newList) <- names
# 			print(names(newList))
# 			setwd(currDir)
# 			return (newList)
# 			}
# 		else {
# 			allLines <-  aggregate(.load.data(dir()[i])$x,  .load.data(dir()[i])["distance"], mean)
# 			newList[[length(newList)+1L]] <-  data.frame(distance = allLines[,1]*40/length(allLines[,1]), x= allLines[,2])
# 			temp <- paste(substr(basename(dir()[i]),1,numDig), "", sep="")
# 			names[i] <- strsplit(temp,".txt")[[1]][1]
# 			getData(i+1, newList, names)
# 		}
# 	}
# 	setwd(directoryPath)
# 	i <-1
# 	names <- c()
# 	findMin <- c()
# 	getData(i, newList, names)
# }
