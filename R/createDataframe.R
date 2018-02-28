#' Dataframe creation

#' @description \code{createDataframe} saves the calculated resistance, perseverence and sensitivity estimates

#' @inheritParams maxLik
#' @param nameVector either a logial value or a character vector. Supported values are \code{nameVector} = "TRUE" to assign the photograph name to the 'name' column, \code{nameVector} = "FALSE" to assign th photograph number to the 'name' column, or \code{nameVector} = a vector the same length as the number of photographs indicating the desired names.
#' @param typeVector a logical value. \code{typeVector} = "TRUE" will add a 'type' vector to the dataframe using values found in the \code{typePlace} position of the photograph names (see \code{\link{IJMacro}} for more details) while \code{typeVector} = "FALSE" will not add a type column.
#' @param typePlace a number that indicates the position of the photograph name to be stored as the 'type' vector'. Defaults to 2. For more details see \code{\link{IJMacro}}
#' @param typeName a character string that indicates what to name the typeVector. Defaults to "type".
#' @param removeClear a logical value that indicates whether to remove the clear halo picture from the dataset (i.e., is this picture an experimental picture, or one solely included to use as a clear halo). Defaults to FALSE.
#' @param standType either "one" or "indiv" to determine whether to use one standard for all photos or individually standardize each photo. Note that "indiv" standardizations are not compatible with measuring FoG.
#' @param needMap Is there a coordinates map to use to assign drug names. This is used when the plates have 16 disks and \code{\link{IJMacro16}} was used rather than \code{\link{IJMacro}. Defaults to "FALSE".
#' @param addZOI Automatically calculate the ZOI from RAD values (RAD*2). Defaults to "TRUE".
#' @param needZOI For standType "indiv", indicates whether to calculate. Defaults to "FALSE". For standType "one" FoG is automatically calculated.
#'@param RADcrit What is the critical inhibition point to use to calculate the ZOI? Defaults to 20% (RAD80), other acceptable values are 50% (RAD50) or 80% (RAD20).
#'@param addSIR Will determine whether the ZOI values match CLSI 'susceptible', 'intermediate', or 'resistant' values.  Currently only implemented for plates with 16 disks (i.e., when \code{\link{IJMacro16}} was run.

#' @details A dataframe with 11 columns:
#' \itemize{
#' 		\item\bold{name:} determined by \code{nameVector}, either photograph names, photograph numbers, or a user-supplied list of names
#'	 	\item\bold{line:} the first components of the \code{namesVector}; everything that comes before the first "_" in the photograph name
#' 		\item\bold{type:} the location within the \code{name} of the phograph type is supplied by \code{typePlace}. Use \code{\link{addType}} if more than one type column are desired.
#' 		\item\bold{RAD80, RAD50, RAD20:} resistance parameters, coresponding to the distance in mm of 80\%, 50\% and 20\% reduction in growth
#' 		\item\bold{FoG80, FoG50, FoG20:} perseverence parameters, coresponding to the fraction of growth achieved above the 80\%, 50\% and 20\% reduction in growth points
#' 		\item\bold{slope:} sensitivity parameter, indicating the slope at the midpoint, i.e., how rapidly the population changes from low growth to full growth
#'	}

#' @return A dataframe "projectName.df" is saved to the global environment and a .csv file "projectName_df.csv" is exported to the "parameter_files" directory.

#' @examples
#' \dontrun{
#' createDataframe("myProject", clearHalo=1)
#' createDataframe("myProject", clearHalo=1, removeClear = TRUE, typeName = "drugAmt")
#' }

#' @export


createDataframe <- function(projectName, clearHalo, diskDiam = 6, maxDist = 30, RADcrit = "20%", standardLoc = 2.5, removeClear = FALSE, nameVector=TRUE, typeVector=TRUE, typePlace=2, typeName = "type", needMap = FALSE, standType = "one", addZOI = TRUE, addSIR=FALSE,  needFoG=FALSE, RADcrit = "20%"){
if(standType=="one"){
	if(!(hasArg(clearHalo))){
		cont <- readline(paste("Please specify photograph number with a clear halo ", sep=""))
		clearHalo <- as.numeric(cont)
		}
	}
	data <- eval(parse(text=projectName))
	if(needMap){
		 mapDir <- file.path(getwd(), "disk_coordinates", projectName)
		 map <- read.csv(file.path(mapDir, paste0(projectName, "_ResultsTable.txt")), sep="\t")
		 photoNames <- unique(unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1])))
		 drugPos <- c()
			for(m in photoNames){
				temp <- subset(map, photoName == m)
				drugPos <- append(drugPos, temp$drug[as.numeric(sort(as.character(temp$XYpos)))])
			}
 		}

	df <- data.frame()
	dotedge <- diskDiam/2 + 0.7
	newdir <- file.path(getwd(), "parameter_files")
	newdir2 <- file.path(getwd(), "parameter_files", projectName)
	newdir3 <- file.path(getwd(), "figures", projectName)

	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df.csv", sep=""))

	if (!file.exists(newdir)){
		dir.create(newdir, showWarnings = FALSE)
		cat(paste("\n\tCreating new directory: ", newdir), sep="")
		}
	if (!file.exists(newdir2)){
		dir.create(newdir2, showWarnings = FALSE)
		cat(paste("\nCreating new directory: ", newdir2), sep="")
		}
	if (!file.exists(newdir3)){
		dir.create(newdir3, showWarnings = FALSE)
		cat(paste("\nCreating new directory: ", newdir3), sep="")
		}
	df <- data.frame(row.names = seq(1, length(data)))

	ML <- paste(projectName, ".ML", sep="")
	ML2 <- paste(projectName, ".ML2", sep="")
	ML <- eval(parse(text=ML))
	ML2 <- eval(parse(text=ML2))

	if(standType == "one"){
		dotMax <- max(sapply(data, function(x) {x[which(x[,1] > standardLoc)[1], 2]}))
		stand <-c( sapply(data, function(x) {dotMax-x[which(x[,1] > standardLoc)[1], 2]}))
		clearHaloData <- data[[clearHalo]]
		startX <- which(clearHaloData[,1] > dotedge+0.5)[1]
		stopX <- which(clearHaloData[,1] > maxDist - 0.5)[1]
		clearHaloData <- clearHaloData[startX:stopX, 1:2]
		clearHaloData$x <- clearHaloData$x + stand[clearHalo]
		clearHaloData$distance <- clearHaloData$distance - (dotedge+0.5)
		clearHaloStand <- clearHaloData[1,2]

		slope <- sapply(c(1:length(data)), .findSlope, data=data, ML=ML, ML2 = ML2, stand = stand, dotedge = dotedge, maxDist = maxDist, clearHaloStand = clearHaloStand, standType = "one")

		FoG.df <-  sapply(c(1:length(data)), .findFoG, data=data, ML=ML, ML2 = ML2, stand = stand, dotedge = dotedge,  maxDist = maxDist, clearHaloStand = clearHaloStand, standardLoc = standardLoc)

		x80 <- unlist(FoG.df[1,])
		x50 <- unlist(FoG.df[2,])
		x20 <- unlist(FoG.df[3,])
		FoG80 <- unlist(FoG.df[4,])
		FoG50 <- unlist(FoG.df[5,])
		FoG20 <- unlist(FoG.df[6,])
		maxFoG <- unlist(FoG.df[7,])
		maxFoG80 <- unlist(FoG.df[8,])
		maxFoG50 <- unlist(FoG.df[9,])
		maxFoG20 <- unlist(FoG.df[10,])

		FoG80[slope < 5] <- NA
		FoG50[slope < 5] <- NA
		FoG20[slope < 5] <- NA
		x80[slope < 5] <- 1
		x50[slope < 5] <- 1
		x20[slope < 5] <- 1

		aveFoG80 <- FoG80/x80
		aveFoG50 <- FoG50/x50
		aveFoG20 <- FoG20/x20

	param <- data.frame(RAD80 =round(x80, digits=3), RAD50 = round(x50, digits=3), RAD20 = round(x20, digits=3), FoG80 = round(FoG80/maxFoG80, digits=2), FoG50 = round(FoG50/maxFoG50, digits=2), FoG20 = round(FoG20/maxFoG20, digits=2), slope=round(slope, digits=1))
}

if(standType == "indiv"){

  RAD.df <-  sapply(c(1:length(data)), .findRAD, data=data, ML=ML, ML2 = ML2, dotedge = dotedge,  maxDist = maxDist)
	x80 <- unlist(RAD.df[1,])
	x50 <- unlist(RAD.df[2,])
	x20 <- unlist(RAD.df[3,])
	asym <- unlist(RAD.df[4,])

  slopes <- sapply(1:length(data), .findSlope, data=data, ML=ML, ML2 = ML2, stand = stand, dotedge = dotedge, maxDist = maxDist, standType = "indiv")
	 param <- data.frame(RAD80 = x80, RAD50 = x50, RAD20 = x20, slope = round(unlist(slopes), digits=1))


	if(needFoG){
	  	FoG.df <-  sapply(c(1:length(data)), .findFogIndiv, data=data, ML=ML, ML2 = ML2,  dotedge = dotedge,  maxDist = maxDist)
  	  x80 <- unlist(FoG.df[1,])
  		x50 <- unlist(FoG.df[2,])
  		x20 <- unlist(FoG.df[3,])
  		FoG80 <- unlist(FoG.df[4,])
  		FoG50 <- unlist(FoG.df[5,])
  		FoG20 <- unlist(FoG.df[6,])
  		maxFoG <- unlist(FoG.df[7,])
  		maxFoG80 <- unlist(FoG.df[8,])
  		maxFoG50 <- unlist(FoG.df[9,])
  		maxFoG20 <- unlist(FoG.df[10,])

	  	param <- data.frame(RAD80 =round(x80, digits=3), RAD50 = round(x50, digits=3), RAD20 = round(x20, digits=3), FoG80 = round(FoG80/maxFoG80, digits=2), FoG50 = round(FoG50/maxFoG50, digits=2), FoG20 = round(FoG20/maxFoG20, digits=2), slope=round(unlist(slopes), digits=1))
	}
}

if(needMap){
	label <- paste(names(data), drugPos, sep="-")
	df <- data.frame(name = names(data), photo=unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1])), dfPos = 1:length(data), photoPos = unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][2])), drug = drugPos, df, param)
	df <- df[order(df$photo, as.numeric(df$drug)),]
	}
else{
	if (is.logical(nameVector)){
		if (nameVector){
			line <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1]))
			df <- data.frame(name = names(data), line)
			}
		if (!nameVector){
			line <- seq(1, length(data))
			df <- data.frame(name = names(data), line, df)
		}
	}
	if (!is.logical(nameVector)){
				line <- nameVector
				names <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1]))
				df <- data.frame(names=names, line=line, df)
			}
	if (typeVector){
				type <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][typePlace]))
				df <- data.frame(df, type, param)
			}
	else {
			df$type <- 1
			df <- data.frame(df, param)
		}
	names(df)[3] <- typeName
}

	if(standType == "one"){
		df <- df[order(df$line),]
		df$FoG80[df$FoG80 >1] <- 1
		df$FoG50[df$FoG50 >1] <- 1
		df$FoG20[df$FoG20 >1] <- 1
		df$FoG80[df$RAD80 == 0] <- NA
		df$FoG50[df$RAD50 == 0] <- NA
		df$FoG20[df$RAD20 == 0] <- NA
		df$FoG80[df$RAD80 == 1] <- NA
		df$FoG50[df$RAD50 == 1] <- NA
		df$FoG20[df$RAD20 == 1] <- NA

		if (removeClear)	df <- df[-clearHalo,]
}

	if(addZOI){
		if(RADcrit = "20%")
			df$ZOI <- round(df$RAD80*2+diskDiam, 0)
			df$ZOI[df$RAD80 ==0] <- diskDiam
		}
		if(RADcrit = "50%")
			df$ZOI <- round(df$RAD50*2+diskDiam, 0)
			df$ZOI[df$RAD80 ==0] <- diskDiam
		}
		if(RADcrit = "80%")
			df$ZOI <- round(df$RAD20*2+diskDiam, 0)
			df$ZOI[df$RAD80 ==0] <- diskDiam
		}
		df$ZOI[df$slope < 0 ] <- diskDiam
	}

if(addSIR){
	if(.Platform$OS.type=="windows"){
	  df <- df[order(df$photo, df$drug),]
	  drugFile <- file.path(.libPaths(), "diskImageR", "drugCutoffs.csv")[1]
		drugFile <- gsub("Program Files", "progra~1", drugFile)
		drugCutoffs <- read.csv(drugFile)
	}
	else{
	  df <- df[order(df$photo, df$drug),]
		drugCutoffs <- read.csv(file.path(.libPaths(), "diskImageR", "drugCutoffs.csv"))
	}
	for(i in 1:length(df$drug)){
		if(df$drug[i] %in% drugCutoffs$Drug_abbrev){
			if(df$ZOI[i] <= subset(drugCutoffs, Drug_abbrev == df$drug[i])$Resistant){
				 df$category[i] <- "R"
			 }
			else{
				if(df$ZOI[i] >= subset(drugCutoffs, Drug_abbrev == df$drug[i])$Susceptible) df$category[i] <- "S"
				else df$category[i] <- "I"
			}
		}
	else df$category[i] <- NA
	}
}

	write.csv(df, file=filename, row.names=FALSE)

	dfName <- paste(projectName, ".df", sep="")
	cat(paste("\n", dfName, " has been written to the global environment", sep=""))
	cat(paste("\nSaving file: ", filename,  sep=""))
	cat(paste("\n", projectName, "_df.csv can be opened in MS Excel.",  sep=""))
	assign(dfName, df, inherits=TRUE)
	}

.findFogIndiv <- function(data, ML, ML2, dotedge = 3.4, maxDist = 35, i){
  startX <- which(data[[i]][,1] > dotedge)[1]
	stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
  minD <- min(data[[i]][startX:stopX, "x"])
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  data[[i]]$x <- data[[i]]$x -min(data[[i]]$x)

  xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  yy<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx)

	xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
	yy <- (yy+min(data[[i]]$x))
	yy[yy < 0] <- 0

	x80 <- xx[which.max(yy> max(yy) * 0.8)]
	x50 <- xx[which.max(yy> max(yy) * 0.5)]
	x20 <- xx[which.max(yy> max(yy) * 0.2)]

  xx80 <- xx[xx<x80]
	yy80 <- yy[xx<x80]
	if(length(xx80) == 0){
	  xx80 <- xx[1:2]
	  yy80 <- yy[1:2]
	}
	xx50 <- xx[xx<x50]
	yy50 <- yy[xx<x50]
	if(length(xx50) == 0){
	  xx50 <- xx[1:2]
	  yy50 <- yy[1:2]
	}
	xx20 <- xx[xx<x20]
	yy20 <- yy[xx<x20]
	if(length(xx20) == 0){
	  xx20 <- xx[1:2]
	  yy20 <- yy[1:2]
	}

	id <- order(xx)
	id80 <- order(xx80)
	id50 <- order(xx50)
	id20 <- order(xx20)

		maxFoG <- sum(diff(xx[id])*zoo::rollmean(yy[id], 2))
		maxFoG80 <- exp(x80)*max(yy80)-min(exp(xx80))*max(yy80)
		maxFoG50 <- exp(x50)*max(yy50)-min(exp(xx50))*max(yy50)
		maxFoG20 <- exp(x20)*max(yy20)-min(exp(xx20))*max(yy20)

		FoG80 <- sum(diff(exp(xx80[id80]))*zoo::rollmean(yy80[id80], 2))
		FoG50 <- sum(diff(exp(xx50[id50]))*zoo::rollmean(yy50[id50], 2))
		FoG20 <- sum(diff(exp(xx20[id20]))*zoo::rollmean(yy20[id20], 2))

		 param <- data.frame(x80 = round(exp(x80), digits=3), x50 = round(exp(x50), digits=3), x20 = round(exp(x20), digits=3) , FoG80 = round(FoG80, digits=0), FoG50= round(FoG50, digits=0), FoG20= round(FoG20, digits=0), maxFoG = round(maxFoG, digits=0), maxFoG80 = round(maxFoG80, digits=0), maxFoG50 = round(maxFoG50, digits=0), maxFoG20 = round(maxFoG20, digits=0))

		 # if (exp(param$x80)<1) 	param$x80 <- 0
		 # if (exp(param$x50)<1)	param$x50 <- 0
		 # if (exp(param$x20)<1)	param$x20 <- 0
		 return(param)
		}

.findFoG <- function(data, ML, ML2, stand, clearHaloStand, dotedge = 3.4, maxDist = 35, standardLoc = 2.5, i){
		startX <- which(data[[i]][,1] > dotedge+0.5)[1]
		stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
		data[[i]] <- data[[i]][startX:stopX, 1:2]
		data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand
		data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)

		xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
		yy<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx)
		ploty <- data[[i]]$x
		ploty[ploty < 0] <-0
		asym <- (ML[[i]]$par[1]+min(data[[i]]$x))

		xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
		yy <- (yy+min(data[[i]]$x))
		yy[yy < 0] <- 0
	   if(max(yy) < asym*0.8){
  		x80 <- xx[which.max(yy> max(yy) * 0.8)]
  		x50 <- xx[which.max(yy> max(yy) * 0.5)]
  		x20 <- xx[which.max(yy> max(yy) * 0.2)]
		 }

 		else{
   		x80 <- xx[which.max(yy> asym * 0.8)]
   		x50 <- xx[which.max(yy> asym * 0.5)]
   		x20 <- xx[which.max(yy> asym * 0.2)]
 		}
  	if (x80 < x50) x80 <- xx[which.max(yy> yy[length(yy)] * 0.8)]

		if(exp(x80)>1) xx80 <- seq(log(data[[i]]$distance[1]), log(round(exp(x80))), length=200)
		else xx80 <- seq(log(data[[i]]$distance[1]), log(data[[i]]$distance[2]), length=200)

		if(exp(x50)>1) xx50 <- seq(log(data[[i]]$distance[1]), log(round(exp(x50))), length=200)
		else xx50 <- seq(log(data[[i]]$distance[1]), log(data[[i]]$distance[2]), length=200)

		if(exp(x20)>1) xx20 <- seq(log(data[[i]]$distance[1]), log(round(exp(x20))), length=200)
		else xx20 <- seq(log(data[[i]]$distance[1]), log(data[[i]]$distance[2]), length=200)

		yy <- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx)
		yy80 <- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx80)
		yy50<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx50)
		yy20<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx20)

		yy <- (yy+min(data[[i]]$x))
		yy[yy < 0] <- 0.1
		yy80 <- (yy80+min(data[[i]]$x))
		yy80[yy80 < 0] <- 0.1
		yy50 <- (yy50+min(data[[i]]$x))
		yy50[yy50 < 0] <- 0.1
		yy20 <- (yy20+min(data[[i]]$x))
		yy20[yy20 < 0] <- 0.1

		id <- order(xx)
		id80 <- order(xx80)
		id50 <- order(xx50)
		id20 <- order(xx20)

		maxFoG <- sum(diff(xx[id])*zoo::rollmean(yy[id], 2))
		maxFoG80 <- exp(x80)*max(yy80)
		maxFoG50 <- exp(x50)*max(yy50)
		maxFoG20 <- exp(x20)*max(yy20)

		FoG80 <- sum(diff(exp(xx80[id80]))*zoo::rollmean(yy80[id80], 2))
		FoG50 <- sum(diff(exp(xx50[id50]))*zoo::rollmean(yy50[id50], 2))
		FoG20 <- sum(diff(exp(xx20[id20]))*zoo::rollmean(yy20[id20], 2))

		 param <- data.frame(x80 = round(exp(x80), digits=3), x50 = round(exp(x50), digits=3), x20 = round(exp(x20), digits=3) , FoG80 = round(FoG80, digits=0), FoG50= round(FoG50, digits=0), FoG20= round(FoG20, digits=0), maxFoG = round(maxFoG, digits=0), maxFoG80 = round(maxFoG80, digits=0), maxFoG50 = round(maxFoG50, digits=0), maxFoG20 = round(maxFoG20, digits=0))

		 if (exp(param$x80)<1) 	param$x80 <- 0
		 if (exp(param$x50)<1)	param$x50 <- 0
		 if (exp(param$x20)<1)	param$x20 <- 0
		 return(param)
		}

#Determine the slope
.findSlope <- function(data, ML, ML2, i, stand, clearHaloStand, dotedge = dotedge,  maxDist = maxDist, standType = standType){
  startX <- which(data[[i]][,1] > dotedge)[1]
  stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
	data[[i]] <- data[[i]][startX:stopX, 1:2]
	if(standType == "one") data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand
	if(standType == "indiv") data[[i]]$x <- data[[i]]$x -min(data[[i]]$x[1:20])
	data[[i]]$x[data[[i]]$x < 0] <- 0
	data[[i]]$distance <- data[[i]]$distance - dotedge
	maxY <- min(ML[[i]]$par[1], (ML2[[i]]$par[1]+ML2[[i]]$par[5]))
	allSlope <- summary(lm(data[[i]]$x~data[[i]]$distance))$coefficients[2]
	disk <- which(data[[i]]$x == min(data[[i]]$x[1:20]))[1]
	#maxYplace <- which(data[[i]][disk:length(data[[i]]$x),2] > maxY)[1]+disk
	maxYplace <- which.max(data[[i]][,2])

	if(!is.na(maxYplace[1])){
	   if(maxYplace[1]==1) xxmid <-1:10
	    if(maxYplace[1]!=1){
	      xxmid <- which(data[[i]]$x[disk:length(data[[i]]$x)] > (maxY/2))+disk ###here###
	    }
	  }
	if (is.na(maxYplace[1])) xxmid <- 1:10
	if(length(xxmid)==0){
		      xxmid <- which(data[[i]]$x[disk:length(data[[i]]$x)] > (maxYplace/2))+disk
	}
	if(length(xxmid)==0){
		      xxmid <- 1:10
		    }

  if(xxmid[1] == 1){
		 if(xxmid[10] == 10) midslope <- 10 #changed from [5] == 5
		else midslope <-  xxmid[10]
	}
	if(xxmid[1] != 1) midslope <- xxmid[1]

	if(is.na(maxYplace)){
		 slope <- 0.1
		 if(allSlope > slope) slope <- allSlope
		 return(slope)
	 }
	if(midslope < 10){
		xxSlope <- data[[i]]$distance[1:20]
		yySlope <- data[[i]]$x[1:20]
		slope <- lm(yySlope ~ xxSlope)$coefficients[2]
		if(allSlope > slope) slope <- allSlope
		return(slope)
	}
	if(midslope >= 10){
		xxSlope <- data[[i]]$distance[(midslope-10):(midslope+10)]
		yySlope <- data[[i]]$x[(midslope-10):(midslope+10)]
		yySlope[yySlope<0] <- 0
		slope <- lm(yySlope ~ xxSlope)$coefficients[2]
		if(allSlope > slope) slope <- allSlope
		return(slope)
		}
}

.findRAD <- function(data, ML, ML2, dotedge, maxDist, i){
    startX <- which(data[[i]][,1] > dotedge)[1]
		stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
		data[[i]] <- data[[i]][startX:stopX, 1:2]
		# minD <- min(data[[i]][startX:stopX, "x"])
		data[[i]]$x <- data[[i]]$x -min(data[[i]]$x[1:20])
		data[[i]]$x[data[[i]]$x < 0] <- 0
		data[[i]]$distance <- data[[i]]$distance - dotedge
		asym <- min(ML[[i]]$par[1], (ML2[[i]]$par[1]+ML2[[i]]$par[5]))
		disk <- min(which(data[[i]]$x[1:20] == 0))[1]
		whichX80 <- which(data[[i]]$x > (asym * 0.8))

    if(length(whichX80) == 0) x80 <- data[[i]]$distance[1]
    if(length(whichX80) > 0){
		  if(whichX80[1] != 1) x80 <- data[[i]]$distance[whichX80[1]]
		  if(whichX80[1] == 1) x80 <- data[[i]]$distance[which(data[[i]]$x[disk+1:length(data[[i]][,1])] > asym * 0.8)[1]+disk]
    }

    whichX50 <- which(data[[i]]$x > (asym * 0.5))
    if(length(whichX50) == 0) x50 <- data[[i]]$distance[1]
    if(length(whichX50) > 0){
		  if(whichX50[1] != 1) x50 <- data[[i]]$distance[whichX50[1]]
			if(whichX50[1] == 1) x50 <- data[[i]]$distance[which(data[[i]]$x[disk+1:length(data[[i]][,1])] > asym * 0.5)[1]+disk]
    }

		whichX20 <- which(data[[i]]$x > (asym * 0.2))
		if(length(whichX20) == 0) x20 <- data[[i]]$distance[1]
		if(length(whichX20) > 0){
		  if(whichX20[1] != 1) x20 <- data[[i]]$distance[whichX20[1]]
			if(whichX20[1] == 1) x20 <- data[[i]]$distance[which(data[[i]]$x[disk+1:length(data[[i]][,1])] > asym * 0.2)[1]+disk]
		}

				# if (x80<1 | is.na(x80)) x80 <- 0
				# if (x50<1 | is.na(x50))	x50 <- 0
				# if (x20<1 | is.na(x20))	x20 <- 0

				if ( is.na(x80)) x80 <- 0
				if ( is.na(x50))	x50 <- 0
				if ( is.na(x20))	x20 <- 0

		 # slopeML <- ML[[i]]$par[3]
		 param <- data.frame(x80 = round(x80, digits=3), x50 = round(x50, digits=3), x20 = round(x20, digits=3), asym = round(asym, digits=2))
		 return(param)
		}
