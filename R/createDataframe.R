#' Dataframe creation

#' @description Writes the main dataframe with resistance, tolerance and sensitivity parameters

#' @inheritParams maxLik
#' @param nameVector either a logial value or a character vector. Supported values are \code{nameVector} = "TRUE" to assign the photograph name to the 'name' column, \code{nameVector} = "FALSE" to assign th photograph number to the 'name' column, or \code{nameVector} = a vector the same length as the number of photographs indicating the desired names.
#' @param typeVector a logical value. \code{typeVector} = "TRUE" will add a 'type' vector to the dataframe using values found in the \code{typePlace} position of the photograph names (see \code{\link{IJMacro}} for more details) while \code{typeVector} = "FALSE" will not add a type column.
#' @param typePlace a number that indicates the position of the photograph name to be stored as the 'type' vector'. Defaults to 2. For more details see \code{\link{IJMacro}}
#' @param typeName a character string that indicates what to name the typeVector. Defaults to "type".

#' @details A dataframe with 11 columns:
#' \itemize{
#' 		\item\bold{name:} determined by \code{nameVector}, either photograph names, photograph numbers, or a user-supplied list of names
#'	 	\item\bold{line:} the first components of the \code{namesVector}; everything that comes before the first "_" in the photograph name
#' 		\item\bold{type:} the location within the \code{name} of the phograph type is supplied by \code{typePlace}. Use \code{\link{addType}} if more than one type column are desired.
#' 		\item\bold{ZOI80, ZOI50, ZOI20:} resistance parameters, coresponding to the distance in mm of 80\%, 50\% and 20\% reduction in growth
#' 		\item\bold{fAUC80, fAUC50, fAUC20:} tolerance parameters, coresponding to the fraction of growth achieved above the 80\%, 50\% and 20\% reduction in growth points
#' 		\item\bold{slope:} sensitivity parameter, indicating the slope at the midpoint, i.e., how rapidly the population changes from low growth to full growth
#'	}
	
#' @return A dataframe "projectName.df" is saved to the global environment and a .csv file "projectName_df.csv" is exported to the "parameter_files" directory. 

#' @export

#' @author Aleeza C. Gerstein

#addin removal of blank disk plate
#try to automate clearHalo

createDataframe <- function(projectName, clearHalo, diskDiam = 6, maxDist = 30, nameVector=TRUE, typeVector=TRUE, typePlace=2, typeName = "type"){
	if(!(hasArg(clearHalo))){
		cont <- readline(paste("Please specify photograph number with a clear halo ", sep=""))
		clearHalo <- as.numeric(cont)
	}
	data <- eval(parse(text=projectName))
	df <- data.frame()
	dotedge <- diskDiam/2 + 0.4
	standardLoc <- 2.5
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

	dotMax <- max(sapply(data, function(x) {x[which(x[,1] > standardLoc)[1], 2]})) 		
	stand <-c( sapply(data, function(x) {dotMax-x[which(x[,1] > standardLoc)[1], 2]}))
	clearHaloData <- data[[clearHalo]]
	startX <- which(clearHaloData[,1] > dotedge+0.5)[1]
	stopX <- which(clearHaloData[,1] > maxDist - 0.5)[1]
	clearHaloData <- clearHaloData[startX:stopX, 1:2]
	clearHaloData$x <- clearHaloData$x + stand[clearHalo] 
	clearHaloData$distance <- clearHaloData$distance - (dotedge+0.5)
	clearHaloStand <- clearHaloData[1,2]

	slope <- sapply(c(1:length(data)), .findSlope, data=data, ML=ML, stand = stand, dotedge = dotedge, maxDist = maxDist, clearHaloStand = clearHaloStand)

	AUC.df <-  sapply(c(1:length(data)), .findAUC, data=data, ML=ML, ML2 = ML2, stand = stand, dotedge = dotedge,  maxDist = maxDist, clearHaloStand = clearHaloStand, standardLoc = standardLoc)

	x80 <- unlist(AUC.df[1,])	
	x50 <- unlist(AUC.df[2,])	
	x20 <- unlist(AUC.df[3,])
	AUC80 <- unlist(AUC.df[4,])	
	AUC50 <- unlist(AUC.df[5,])	
	AUC20 <- unlist(AUC.df[6,])	
	maxAUC <- unlist(AUC.df[7,])		
	maxAUC80 <- unlist(AUC.df[8,])		
	maxAUC50 <- unlist(AUC.df[9,])		
	maxAUC20 <- unlist(AUC.df[10,])				
	
	AUC80[slope < 5] <- NA
	AUC50[slope < 5] <- NA
	AUC20[slope < 5] <- NA
	x80[slope < 5] <- 1
	x50[slope < 5] <- 1
	x20[slope < 5] <- 1

	aveAUC80 <- AUC80/x80
	aveAUC50 <- AUC50/x50
	aveAUC20 <- AUC20/x20	

	param <- data.frame(ZOI80 =round(x80, digits=0), ZOI50 = round(x50, digits=0), ZOI20 = round(x20, digits=0), fAUC80 = round(AUC80/maxAUC80, digits=2), fAUC50 = round(AUC50/maxAUC50, digits=2), fAUC20 = round(AUC20/maxAUC20, digits=2), slope=round(slope, digits=1))
	
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

	df <- df[order(df$line),] 
	df$fAUC80[df$fAUC80 >1] <- 1
	df$fAUC50[df$fAUC50 >1] <- 1
	df$fAUC20[df$fAUC20 >1] <- 1	
	df$fAUC80[df$ZOI80 == 1] <- 1
	df$fAUC50[df$ZOI50 == 1] <- 1
	df$fAUC20[df$ZOI20 == 1] <- 1
	
	write.csv(df, file=filename, row.names=FALSE)	
	
	dfName <- paste(projectName, ".df", sep="")
	cat(paste("\n", dfName, " has been written to the global environment", sep=""))
	cat(paste("\nSaving file: ", filename,  sep=""))
	cat(paste("\n", projectName, "_df.csv can be opened in MS Excel.",  sep=""))
	assign(dfName, df, envir=globalenv())
	}

#Determine the slope

.findSlope <- function(data, ML, i, stand, clearHaloStand, dotedge = 3.4,  maxDist = 35){
	startX <- which(data[[i]][,1] > dotedge+0.5)[1]
	stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
	data[[i]] <- data[[i]][startX:stopX, 1:2]
	data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand 
	data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
	xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
	yy<- .curve(ML[[i]]['par'][1]$par[1], ML[[i]]['par'][1]$par[2], ML[[i]]['par'][1]$par[3],xx)
	yycor <- (yy+min(data[[i]]$x))
	xcross <- exp(ML[[i]]['par'][1]$par[2])
	xxmid <- which.max(exp(xx) > xcross)
	if ((xxmid-10) > 1){
		xxSlope <- xx[(xxmid-10):(xxmid+10)]
		yySlope <- yy[(xxmid-10):(xxmid+10)]
		}
	else {
		xxSlope <- xx[1:(xxmid+10)]
		yySlope <- yy[1:(xxmid+10)]
	}
	slope <- lm(yySlope ~ xxSlope)$coefficients[2]
	return(slope)
}

.findAUC <- function(data, ML, ML2, stand, clearHaloStand, dotedge = 3.4, maxDist = 35, standardLoc = 2.5, i){	
	startX <- which(data[[i]][,1] > dotedge+0.5)[1]
	stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
	data[[i]] <- data[[i]][startX:stopX, 1:2]
	data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand 
	data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
	xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200) 
	yy<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx) 
	# ic50 <- ML[[i]]$par[2]	
	ploty <- data[[i]]$x
	ploty[ploty < 0] <-0
	# slope <- ML[[i]]$par[3]
	asym <- (ML[[i]]$par[1]+min(data[[i]]$x))

	xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200) 				
	yy <- (yy+min(data[[i]]$x))
	yy[yy < 0] <- 0		
	x80 <- xx[which.max(yy> asym * 0.2)]
	x50 <- xx[which.max(yy> asym * 0.5)]
	x20 <- xx[which.max(yy> asym * 0.8)]
	if (x20 < x50) x20 <- xx[which.max(yy> yy[length(yy)] * 0.8)]

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

	maxAUC <- sum(diff(xx[id])*zoo::rollmean(yy[id], 2))	
	maxAUC80 <- exp(x80)*max(yy80)
	maxAUC50 <- exp(x50)*max(yy50)	
	maxAUC20 <- exp(x20)*max(yy20)	
		
	AUC80 <- sum(diff(exp(xx80[id80]))*zoo::rollmean(yy80[id80], 2))		
	AUC50 <- sum(diff(exp(xx50[id50]))*zoo::rollmean(yy50[id50], 2))		
	AUC20 <- sum(diff(exp(xx20[id20]))*zoo::rollmean(yy20[id20], 2))		
	
	 param <- data.frame(x80 = round(exp(x80), digits=0), x50 = round(exp(x50), digits=2), x20 = round(exp(x20), digits=0) , AUC80 = round(AUC80, digits=0), AUC50= round(AUC50, digits=0), AUC20= round(AUC20, digits=0), maxAUC = round(maxAUC, digits=0), maxAUC80 = round(maxAUC80, digits=0), maxAUC50 = round(maxAUC50, digits=0), maxAUC20 = round(maxAUC20, digits=0))	
	 
	 if (exp(param$x80)<1){
	 	param$x80 <- 1}
	 if (exp(param$x50)<1){
	 	param$x50 <- 1}	 
	if (exp(param$x20)<1){
		param$x20 <- 1}	  
		return(param)	
	}
