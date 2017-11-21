#' Find the points of inhibition, maximum growth and resource limitations

#' @description \code{maxLik}  uses the intensity information to find the length of inhibition, points of maximum growth, and points of resource limitation.

#' @param projectName the short name in use for the project.
#' @param diskDiam the diameter of the diffusion disk in mm, defaults to 6.
#' @param maxDist a numeric value indicating the maximum distance away from the disk to be considered. Defaults to 30mm.
#' @param xplots a numeric value indicating how many plots to plot in each row, does not influence maximum likelihood fitting
#' @param ymax a numeric value indicating the maximum y value plotted in each graph, does not influence maximum likelihood fitting
#' @param height a numeric value indicating the height of the pdf file generated, does not influence maximum likelihood fitting
#' @param width a numeric value indicating the width of the pdf file generated, does not influence maximum likelihood fitting
#' @param RAD a numeric value indicating the critical level of the radius of inhibition (i.e., resistance) parameter to plot, does not influence maximum likelihood fitting. Currently only \code{RAD} = "80" (80\% reduction in growth), \code{RAD} = "50" (50\% reduction in growth), \code{RAD} = "20" (20\% reduction in growth), and \code{RAD} = "all" are supported.
#' @param popUp a logical value indicating whether to pop up the figure after it has been created.
#' @param nameVector either a logial value indicating whether to plot the photograph names above the graph or not or a vector the same length as the number of pictures containing the desired names. Defaults to TRUE.
#' @param overwrite a logical value indicating whether to overwrite existing figures created on the same day for the same project name.defaults to TRUE.
#' @param plotParam a logical value indicating whether to save plots containing, at minimum, the fitted logistic equatoin and specified RAD levels to plot, but may also include the FoG \code{plotFoG} = "TRUE" or the components of the logistic equation \code{plotCompon} = "TRUE". Defaults to TRUE.
#' @param savePDF a logical value indicating whether to save a PDF file or open a new quartz. Defaults to TRUE.
#' @param plotSub allows you to plot only a subset of photographs - indicate with a vector the corresponding numeric indices of the data you wish to plot. Photographs are numbered alphabetically by name, and the photograph numbers can also be found by using the showNum option in \code{\link{plotRaw}}. Defaults to NA, which will plot data from all photographs. Note this does not affect the analysis component, all data is always analyzed.

#' @details \code{\link{maxLik}} searches for the maximum likelihood (ML) parameter for a single logistic and double logistic equation using the pixel intensity information previously determined from \code{\link{IJMacro}}. The equations fit are
#' single logistic ('ML'):
#'	y = asymA*exp(scalA(x-od50A))\(1+exp(scalA(x-od50A)))+N(0, sigma)
#'  double logistic ('ML2'):
#'	y = asymA*exp(scalA(x-od50A))\(1+exp(scalA(x-od50A)))+asymB*exp(scalB(x-od50B)))\(1+exp(scalB(x-od50B)))+N(0, sigma)
#' where asymA and asymB are the two asymptotes, od50A and odB are the midpoints (of the two curves), scalA and scalB are the slopes at odA and odB divided by asymA/4 and asymB/4, respectively. Specifically, \code{\link{maxLik}} uses the\code{\link[subplex]{subplex}} method of \code{\link{optim}}, as implemented in \code{\link[diversitree]{find.mle}}.  The single logistic is the essentially the same as the double, yet fits only a single asymptote, midpoint and slope. The results from the double logistic fit are used in \code{\link{createDataframe}} to find the resistance points as well as to fit the area under the curve and thus tolerance, the single logistic is used to determine the midpoint of the curve which is later used to find sensitivity, i.e., the slope at this midpoint.

#' @section Warning:
#' Depending on the number of photographs to be analyzed, `maxLik()` can take a fair amount of time, upwards of an hour or more. This is due to the maximum likelihood fitting procedures, which determine the best fit parameters from multiple different starting values. The status is indicated by a series of dots (".") in the R console, with one dot per photograph. If for some reason the procedure gets halted in the middle of \code{maxLik()} (e.g., computer is shut down) as long as R remains open it should resume where it left off when the computer is reactivated.

#' @return Two lists, ML and ML2 are saved to the global environment. A pdf file with one plot for each photograph is saved to visualize the results of curve fitting, zone of inhibition (resistance) and the area under the curve (tolerance).

#' @export

#' @examples
#' \dontrun{
#' maxLik("myProject", clearHalo=1)
#' maxLik("myProject", clearHalo=1, xplots = 2, height = 4, width = 6, needML = FALSE)
#' }

inhibGrowPts <- function(projectName, diskDiam = 12.7, maxDist=30, ymax=125, xplots = 5, height = 8,  width = 8, plotPts="all", popUp = TRUE, nameVector=TRUE, overwrite = TRUE, savePDF= TRUE){
# FoG=20, plotFoG = TRUE,  needML = TRUE,
	options(warn=-1)

	if(!plotPts %in% c(90, 75, 50, 25, "all")){
		stop("Current suppported growth point values for plotting = 'all', 90, 75, 50, 25, 10")
		}

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

	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df.csv", sep=""))
	df <- data.frame(row.names = seq(1, length(data)))

#standard dot edge calculation - takes awhile to get down to dark so diskDiam/2-1
dotedge <- diskDiam/2+1

#where is the minimum intensity - can use for later calculations
minIntensity <- c(sapply(data, function(x) {min(x[min(which(x[,1]>dotedge)):max(which(x[,1]<maxDist)),2])}))

whichMinIntensity <- c(mapply(function(x, y) {which(x[,2]==y)}, x = data, y = minIntensity))
whereMinIntensity <- c(mapply(function(x, y) {x[y, 1]}, x = data, y = whichMinIntensity))
whichMinIntensity5 <- c(mapply(function(x, y) {max(which(x[,2] < minIntensity*1.05))}, x = data, y = minIntensity))
whereMinIntensity5 <- c(mapply(function(x, y) {x[y, 1]}, x = data, y = whichMinIntensity5))

#subtract the minimum intensity from each intensity value
dataCor <- lapply(data, function(x) {x[,2] - min(x[,2])})

#add the "corrected data" vector to the data list
datas <- Map(cbind, data, dataCor)
data <- datas

#where is the maximum intensity? This is also the point of growth maximum
maxIntensity <- c(sapply(data, function(x) {max(x[min(which(x[,1]>dotedge)):max(which(x[,1]<maxDist)), 3])}))
whichMaxIntensity <- mapply(function(x, y) {which(x[,3]==y)}, x = data, y = maxIntensity)
whereMaxIntensity <- mapply(function(x, y) {x[which(x[,3]==y)[1],1]}, x = data, y = maxIntensity)

#what is the slope between minimum intensity (*1.05) and maximum intensity = how sharp is the transition between inhibition and growth?
slope2Max <- round(mapply(function(x, y, z) {coefficients(lm(x[y:z, 2]~ x[y:z, 1]))[2]}, x = data, y = whichMinIntensity5, z = whichMaxIntensity), digits=2)

#use this, change whereMaxIntensity in the param data frame to be corrected for diskDiam
#distance2Max <- whereMaxIntensity-diskDiam/2

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

param <- data.frame(whereMaxIntensity =round(distance2Max, digits=2), maxIntensity = round(maxIntensity, digits=2),  whichMaxIntensity = whichMaxIntensity, whereMinIntensity = round(whereMinIntensity, 2), minIntensity = round(minIntensity, 2), whereMinIntensity5 = whereMinIntensity5, whichMinIntensity5 = whichMinIntensity5,  slope2Max = slope2Max,  dist2_90 = distance2Max+round(dist2_90, digits=2), dist2_75 = distance2Max+round(dist2_75, digits=2), dist2_50 = distance2Max+round(dist2_50, digits=2), dist2_25 = distance2Max+round(dist2_25, digits=2), dist2_10 = distance2Max+round(dist2_10, digits=2), min2_90 = round(min2_90, 2), min2_75 = round(min2_75, 2), min2_50 = round(min2_50, 2), min2_25 = round(min2_25, 2), min2_10 = round(min2_10, 2))

if (is.logical(nameVector)){
	if (nameVector){
		line <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1]))
		df <- data.frame(name = names(data), line, param)
		}

	if (!nameVector){
		line <- seq(1, length(data))
		df <- data.frame(name = names(data), line, df, param)
	}
}
if (!is.logical(nameVector)){
	line <- nameVector
	names <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][1]))
	df <- data.frame(names=names, line=line, df, param)
	}

df <- df[order(df$line),]
write.csv(df, file=filename, row.names=FALSE)

dfName <- paste(projectName, ".df", sep="")
cat(paste("\n", dfName, " has been written to the global environment", sep=""))
cat(paste("\nSaving file: ", filename,  sep=""))
cat(paste("\n", projectName, "_df.csv can be opened in MS Excel.",  sep=""))
# assign(dfName, df, envir=globalenv())
assign(dfName, df, inherits=TRUE)
}
