#' Find the points of inhibition, maximum growth and resource limitations

#' @description \code{inhibGrow}  uses the intensity information to find the length of inhibition, points of maximum growth, and points of resource limitation.

#' @param projectName the short name in use for the project.
#' @param diskDiam the diameter of the diffusion disk in mm, defaults to 6.
#' @param maxDist a numeric value indicating the maximum distance away from the disk to be considered. Defaults to 30mm.
#' @param nameVector either a logial value indicating whether to plot the photograph names above the graph or not or a vector the same length as the number of pictures containing the desired names. Defaults to TRUE.

#' @details \code{\link{inhibGrow}} is useful for compounds that both inhibit growth yet are also utilized. In comparison to \code{\link{maxLik}}, the growth profile for these compounds don't reach an asymptote, but rather are characterized by inhibition, followed by a point of maximum growth and then a growth reduction.
#' The pixel intensity information previously determined from \code{\link{IJMacro}} is used to determine the growth maximum and inhibition points.
#' Slopes are calculated from the minimum intensity to the point of maximum intensity ('slope2Max') and the point of maximum intensity to the edge of the plate ('slopeFromMax').

#' @return A dataframe "projectName_inhib.df" is saved to the global environment and a .csv file "projectName_inhib_df.csv" is exported to the "parameter_files" directory.

#' @export

#' @examples
#' \dontrun{
#' inhibGrow("myProject")
#' maxLik("myProject", diskDiam = 12, maxDist = 25)
#' }

inhibGrow <- function(projectName, diskDiam = 12.7, maxDist=30, nameVector = TRUE){
	data <- eval(parse(text=projectName))
	if (is.logical(nameVector)){
		if (nameVector){label <- names(data)}
		else {label <- rep("", length(data))}
	}
	if (!is.logical(nameVector)) label <- nameVector

	newdir <- file.path(getwd(), "parameter_files")
	newdir2 <- file.path(getwd(), "parameter_files", projectName)
	newdir3 <- file.path(getwd(), "figures", projectName)

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

	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_inhib_df.csv", sep=""))
	df <- data.frame(row.names = seq(1, length(data)))

  #standard dot edge calculation - takes awhile to get down to dark so diskDiam/2-1
  dotedge <- diskDiam/2+0.7

  #subtract the minimum intensity from each intensity value
  dataCor <- lapply(data, function(x) {x[,2] - min(x[,2])})
  #add the "corrected data" vector to the data list
  datas <- Map(cbind, data, dataCor)
  data <- datas

  #where is the minimum intensity - used for slope calculation
  minIntensity <- c(sapply(data, function(x) {min(x[1:max(which(x[,1]<maxDist)),2])}))
  whichMinIntensity <- unlist(c(mapply(function(x, y) {which(x[,2]==y)[1]}, x = data, y = minIntensity)))
  whereMinIntensity <- unlist(c(mapply(function(x, y) {x[y, 1]}, x = data, y = whichMinIntensity)))
  whichMinIntensity5 <- c(mapply(function(x, y) {min(which(x[whichMinIntensity:length(x),2] < minIntensity*1.05))}, x = data, y = minIntensity))
  whereMinIntensity5 <- c(mapply(function(x, y) {x[y, 1]}, x = data, y = whichMinIntensity5))+whereMinIntensity

  #where is the maximum intensity? This is also the point of growth maximum
  maxIntensity <- c(sapply(data, function(x){max(x[min(which(x[,1]>dotedge)):max(which(x[,1]<maxDist)), 3])})) 
  whichMaxIntensity <- mapply(function(x, y) {which(x[,3]==y)[1]}, x = data, y = maxIntensity)
  whereMaxIntensity <- mapply(function(x, y) {x[which(x[,3]==y)[1],1]}, x = data, y = maxIntensity)
  whichMaxIntensity95up <- c(mapply(function(x, y) {min(which(x[,3] < maxIntensity*0.95))}, x = data, y = maxIntensity))
  whereMaxIntensity95up <- c(mapply(function(x, y) {x[y, 1]}, x = data, y = whichMaxIntensity95up))
  
  whichMaxIntensity95down <- c(mapply(function(x, y) {min(which(x[whichMaxIntensity:length(whichMaxIntensity),3] < maxIntensity*0.95))}, x = data, y = maxIntensity))+whichMaxIntensity
  whereMaxIntensity95down <- c(mapply(function(x, y) {x[y, 1]}, x = data, y = whichMaxIntensity95down))

  #what is the slope between minimum intensity (*1.05) and maximum intensity = how sharp is the transition between inhibition and growth?
  slope2Max <- round(mapply(function(x, y, z) {coefficients(lm(x[y:z, 2]~ x[y:z, 1]))[2]}, x = data, y = whichMinIntensity, z = whichMaxIntensity), digits=2)
  #slopeFromMax <- round(mapply(function(x, y, z) {coefficients(lm(x[y:which(x[,1]<maxDist, 2])~ x[y:which(x[,1]<maxDist, 1]))[2]}, x = data, y = whichMaxIntensity), digits=2)
  
  #Find the cutoff intensity values and distance for post-max growth points of resource limitation
  min2_90 <- maxIntensity*0.9
  min2_75 <- maxIntensity*0.75
  min2_50 <- maxIntensity*0.5
  min2_25 <- maxIntensity*0.25
  min2_10 <- maxIntensity*0.1

  #note that these are all distances past the point of maximal growth
  dist2_90 <- c(mapply(function(x, y, z) {x[min(which(x[which(x[,3]==y):max(which(x[,1]<maxDist)), 3] < z)),1]}, x = data, y = maxIntensity, z = min2_90))
  dist2_75 <- c(mapply(function(x, y, z) {x[min(which(x[which(x[,3]==y):max(which(x[,1]<maxDist)), 3] < z)),1]}, x = data, y = maxIntensity, z = min2_75))
  dist2_50 <- c(mapply(function(x, y, z) {x[min(which(x[which(x[,3]==y):max(which(x[,1]<maxDist)), 3] < z)),1]}, x = data, y = maxIntensity, z = min2_50))
  dist2_25 <- c(mapply(function(x, y, z) {x[min(which(x[which(x[,3]==y):max(which(x[,1]<maxDist)), 3] < z)),1]}, x = data, y = maxIntensity, z = min2_25))
  dist2_10 <- c(mapply(function(x, y, z) {x[min(which(x[which(x[,3]==y):max(which(x[,1]<maxDist)), 3] < z)),1]}, x = data, y = maxIntensity, z = min2_10))


  distance2Max <- whereMaxIntensity-diskDiam/2

  #legacy with minIntensity values
  #param <- data.frame(whereMaxIntensity =round(distance2Max, digits=2), maxIntensity = round(maxIntensity, digits=2),  whichMaxIntensity = unlist(whichMaxIntensity), whereMinIntensity = round(whereMinIntensity, 2)-diskDiam/2, minIntensity = round(minIntensity, 2), whereMinIntensity5 = whereMinIntensity5-diskDiam/2, whichMinIntensity5 = whichMinIntensity5,  slope2Max = slope2Max,  dist2_90 = distance2Max+round(dist2_90, digits=2), dist2_75 = distance2Max+round(dist2_75, digits=2), dist2_50 = distance2Max+round(dist2_50, digits=2), dist2_25 = distance2Max+round(dist2_25, digits=2), dist2_10 = distance2Max+round(dist2_10, digits=2), min2_90 = round(min2_90, 2), min2_75 = round(min2_75, 2), min2_50 = round(min2_50, 2), min2_25 = round(min2_25, 2), min2_10 = round(min2_10, 2))
#slopeFromMax = slopeFromMax,
  paramKeep <- data.frame(whereMaxIntensity =round(distance2Max, digits=2), maxIntensity = round(maxIntensity, digits=2), slope2Max = slope2Max,  dist2_90 = distance2Max+round(dist2_90, digits=2), dist2_75 = distance2Max+round(dist2_75, digits=2), dist2_50 = distance2Max+round(dist2_50, digits=2), dist2_25 = distance2Max+round(dist2_25, digits=2), dist2_10 = distance2Max+round(dist2_10, digits=2), min2_90 = round(min2_90, 2), min2_75 = round(min2_75, 2), min2_50 = round(min2_50, 2), min2_25 = round(min2_25, 2), min2_10 = round(min2_10, 2))

  if (is.logical(nameVector)){
  	if (nameVector){
  		line <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1]))
  		df <- data.frame(name = names(data), line, paramKeep)
  	}
	  if (!nameVector){
  		line <- seq(1, length(data))
  		df <- data.frame(name = names(data), line, df, paramKeep)
	  }
  }
  if (!is.logical(nameVector)){
  	line <- nameVector
  	names <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1]))
  	df <- data.frame(names=names, line=line, df, paramKeep)
	}

  df <- df[order(df$line),]
  write.csv(df, file=filename, row.names=FALSE)

  dfName <- paste(projectName, "_inhib.df", sep="")
  cat(paste("\n", dfName, " has been written to the global environment", sep=""))
  cat(paste("\nSaving file: ", filename,  sep=""))
  cat(paste("\n", projectName, "_inhib_df.csv can be opened in MS Excel.",  sep=""))
  # assign(dfName, df, envir=globalenv())
  assign(dfName, df, inherits=TRUE)
}