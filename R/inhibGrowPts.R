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
minIntensity <- c(sapply(data, function(x) {min(x[,2])}))
whichMinIntensity <- c(sapply(data, function(x) {which(x[,2]==min(x[,2]))}))
whereMinIntensity <- c(sapply(data, function(x) {x[whichMinIntensity, 1]}))
whichMinIntensity5 <- c(sapply(data, function(x) {max(which(x[,2]<min(x[,2]*1.05)))}))
whereMinIntensity5 <- c(sapply(data, function(x) {x[whichMinIntensity5, 1]}))

#subtract the minimum intensity from each intensity value
dataCor <- lapply(data, function(x) {x[,2] - min(x[,2])})

#add the "corrected data" vector to the data list
datas <- Map(cbind, data, dataCor)
data <- datas

#where is the maximum intensity? This is also the point of growth maximum
maxIntensity <- c(sapply(data, function(x) {max(x[min(which(x[,1]>dotedge)):max(which(x[,1]<maxDist)), 3])}))
whichMaxIntensity <- mapply(function(x, y) {which(x[,3]==y)}, x = data, y = maxIntensity)
whereMaxIntensity <- round(c(mapply(function(x, y) {x[which(x[,3]==y),1]}, x = data, y = maxIntensity)), digits=2)

#what is the slope between minimum intensity (*1.05) and maximum intensity = how sharp is the transition between inhibition and growth?
slope2Max <- round(mapply(function(x, y, z) {coefficients(lm(x[y:z, 2]~ x[y:z, 1]))[2]}, x = data, y = whichMinIntensity5, z = whichMaxIntensity), digits=2)

distance2Max <- whereMaxIntensity-diskDiam/2

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

param <- data.frame(GRMax.X =round(distance2Max, digits=2), GRMax.Y = round(maxIntensity, digits=2), slope2Max = slope2Max,  dist2_90 = distance2Max+round(dist2_90, digits=2), dist2_75 = distance2Max+round(dist2_75, digits=2), dist2_50 = distance2Max+round(dist2_50, digits=2), dist2_25 = distance2Max+round(dist2_25, digits=2), dist2_10 = distance2Max+round(dist2_10, digits=2))

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

plotParam <- function(x){
	plot(data[[x]][,1], data[[x]][,2], xlab="", xaxt="n", yaxt="n", ylim=c(0, 260))
	arrows(whereMaxIntensity[x], maxIntensity[x]+minIntensity[x]+80, whereMaxIntensity[x], maxIntensity[x]+minIntensity[x]+20, length=0.1)
  text(whereMaxIntensity[x], maxIntensity[x]+minIntensity[x]+100, "maxGR", cex=0.6)
  arrows(whereMaxIntensity[x]+dist2_90[x], min2_90[x]+minIntensity[x]+60, whereMaxIntensity[x]+dist2_90[x], min2_90[x]+minIntensity[x]+20, length=0.05, col="grey")
 arrows(whereMaxIntensity[x]+dist2_75[x], min2_75[x]+minIntensity[x]+60, whereMaxIntensity[x]+dist2_75[x], min2_75[x]+minIntensity[x]+20, length=0.05, col="grey")
 arrows(whereMaxIntensity[x]+dist2_50[x], min2_50[x]+minIntensity[x]+60, whereMaxIntensity[x]+dist2_50[x], min2_50[x]+minIntensity[x]+20, length=0.05, col="grey")
 arrows(whereMaxIntensity[x]+dist2_25[x], min2_50[x]+minIntensity[x]+60, whereMaxIntensity[x]+dist2_25[x], min2_25[x]+minIntensity[x]+20, length=0.05, col="grey")
 # text(whereMaxIntensity[x]+dist2_90[x]-2, min2_90[x]+minIntensity[x]+90, "90%", cex=0.5, pos=4)
 # text(whereMaxIntensity[x]+dist2_75[x]-2, min2_75[x]+minIntensity[x]+70, "75%", cex=0.5, pos=4)
 # text(whereMaxIntensity[x]+dist2_50[x]-2, min2_50[x]+minIntensity[x]+50, "50%", cex=0.5, pos=4)

	segments(whereMinIntensity5[x], lm(data[[x]][whichMinIntensity5[x]:whichMaxIntensity[x], 2]~ data[[x]][whichMinIntensity5[x]:whichMaxIntensity[x], 1])$fitted.values[1], whereMaxIntensity[x], lm(data[[x]][whichMinIntensity5[x]:whichMaxIntensity[x], 2]~ data[[x]][whichMinIntensity5[x]:whichMaxIntensity[x], 1])$fitted.values[length(lm(data[[x]][whichMinIntensity5[x]:whichMaxIntensity[x], 2]~ data[[x]][whichMinIntensity5[x]:whichMaxIntensity[x], 1])$fitted.values)], col="red", lwd=2)

}

par(mfrow=c(2, 1), mar=c(3, 3, 1, 1))
plotParam(1)
axis(1, labels=FALSE)
axis(2, las=2)

plotParam(2)
axis(1)
axis(2, las=2)

#plotting to show it's working
# plot(data[[1]][,1], data[[1]][,3])
# abline(v=whereMax[1])
# abline(v=whereMax[1]+dist2_90[1])
# abline(v=whereMax[1]+dist2_75[1])
# abline(v=whereMax[1]+dist2_50[1])
# abline(v=whereMax[1]+dist2_25[1])
# abline(v=whereMax[1]+dist2_10[1])
#
# plot(data[[2]][,1], data[[2]][,3])
# abline(v=whereMax[2])
# abline(v=whereMax[2]+dist2_90[2])
# abline(v=whereMax[2]+dist2_75[2])
# abline(v=whereMax[2]+dist2_50[2])
# abline(v=whereMax[2]+dist2_25[2])
# abline(v=whereMax[2]+dist2_10[2])

#need to save data back to working space
#plots with points on them
#save a dataframe with the values in it - use this function to do all the extensions for this pupose


####STOPPED EDITING HERE
if(plotParam){
		clearHaloData <- data[[clearHalo]]
		startX <- which(clearHaloData[,1] > dotedge+0.5)[1]
		stopX <- which(clearHaloData[,1] > maxDist - 0.5)[1]
		clearHaloData <- clearHaloData[startX:stopX, 1:2]
		clearHaloData$x <- clearHaloData$x + stand[clearHalo]
		clearHaloData$distance <- clearHaloData$distance - (dotedge+0.5)
		clearHaloStand <- clearHaloData[1,2]

.plotParam(projectName, MLtip=MLtip, dotedge = dotedge, stand = stand, standardLoc = standardLoc, maxDist = maxDist, ymax = ymax, clearHaloStand = clearHaloStand, RAD=RAD, height = height, width=width, xplots = xplots,label=label, overwrite = overwrite, popUp = popUp, savePDF = savePDF, plotSub = plotSub, plotCompon=plotCompon)

		# .plotParam(projectName, ML=ML, ML2=ML2, dotedge = dotedge, stand = stand, standardLoc = standardLoc, maxDist = maxDist, ymax = ymax, clearHaloStand = clearHaloStand, FoG=FoG, RAD=RAD, height = height, width=width, xplots = xplots,label=label, overwrite = overwrite, popUp = popUp, plotFoG = plotFoG, savePDF = savePDF, plotSub = plotSub, plotCompon=plotCompon)
	}
	alarm()
}

.curve <-  function(asym, ic50,scal, x) {asym*exp(scal*(x-ic50))/(1+exp(scal*(x-ic50)))}

.curve2 <- function(asym, od50, scal, asymB, od50B, scalB, x) { asym*exp(scal*(x-od50))/(1+exp(scal*(x-od50)))+asymB*exp(scalB*(x-od50B))/(1+exp(scalB*(x-od50B)))}


.plotParam <- function(projectName, ML , ML2, stand,  clearHaloStand, standardLoc = 2.5, ymax=200, dotedge = 3.4, maxDist= 40, xplots = 4, height = 10, width=7,  FoG=50, RAD=50, overwrite = TRUE, popUp = TRUE, plotFoG = TRUE, label=label, savePDF = TRUE, plotSub = plotSub, plotCompon=plotCompon){
	data <- eval(parse(text=projectName))
	if(is.na(plotSub[1])){
		plotSub <- 1:length(data)
		}
	fileFolder <- projectName
	dir.create(file.path(getwd(), "figures"), showWarnings= FALSE)
	dir.create(file.path(getwd(), "figures", fileFolder), showWarnings= FALSE)
	t <- file.path("figures", projectName , paste(projectName, "_FoG.pdf", sep=""))
	if (!overwrite){
		if (file.exists(t)){
			t <- file.path("figures", projectName , paste(projectName, "_FoG_2_FoG", FoG, "_RAD", RAD, ".pdf", sep=""))
			if (file.exists(t)){
				k <- 2
				while(file.exists(t)){
					k <- k+1
					t <- file.path("figures", projectName, paste(projectName, "_FoG_", k, "_FoG", FoG, "_RAD", RAD, ".pdf", sep=""))
					}
				}
			}
		}

	if(xplots > length(plotSub)){
		xplots <- length(plotSub)
	}
	if (ceiling(length(plotSub)/xplots) < 6) {
		yplots<- ceiling(length(plotSub)/xplots)}
	else {yplots<- 6}
	numpages <- ceiling(length(plotSub)/(xplots*yplots))
	if(savePDF){
		pdf(t, width=width, height=height)
	}
	# if(!savePDF){
		# quartz(width=width, height=height)
	# }
	par(mfrow=c(yplots , xplots), mar=c(1,1,1,1), oma=c(4,5,1,1))
	for (k in plotSub){
			.singleFoG(data = data, ML = ML, ML2 = ML2, dotedge = dotedge, maxDist = maxDist, ymax = ymax, stand = stand, i = k,FoG=FoG, RAD = RAD, clearHaloStand = clearHaloStand, label=label[k], plotFoG = plotFoG, plotCompon=plotCompon)
		if(numpages == 1){
			if (k >= xplots*yplots-xplots+1){
				axis(1, cex.axis=1)
				}
			else {axis(1, cex.axis=1, labels= FALSE)}
			}
		if(numpages == 2){
			if (k >= xplots*yplots-xplots+1 & k < xplots*yplots+1){
				axis(1, cex.axis=1)
				}
			if (k >= 2*xplots*yplots-xplots+1){
				axis(1, cex.axis=1)
				}
			else {axis(1, cex.axis=1, labels= FALSE)}
			}
		if(numpages == 3){
			if (k >= xplots*yplots-xplots+1 & k < xplots*yplots+1 | k >= 2*xplots*yplots-xplots+1 & k < 2*xplots*yplots+1 | k >= 3*xplots*yplots-xplots+1){
				axis(1, cex.axis=1)
				}
			else{axis(1, labels=FALSE)}
			}
		axis(1, labels=FALSE)
		j <- 1
		while (j <= numpages){
			if (k %in% seq(1, j*yplots*xplots, by=xplots)) {axis(2, cex.axis=1, las=2)}
			j <- j+1
		}
	}

	mtext("Distance (mm)", outer=TRUE, side=1, line=2, cex=1.2)
	mtext("Pixel intensity", outer=TRUE, side=2, line=2, cex=1.2)

	if(savePDF){
		dev.off()
		cat(paste("\nFigure saved: ", t, sep=""))

		if(popUp){
			tt <- paste("open", t)
			system(tt)
		}
	}
}
