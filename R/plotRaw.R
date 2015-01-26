#' Used to plot the results of the imageJ analysis

#' @description This function creates a pdf figure of plots showing the results of the imageJ analysis, with one plot for each picture. This function is primarily for visualization purposes.

#' @param projectName the short name to be used for the project
#' @param ymin a numeric value indicating the minimum y value plotted in each graph
#' @param ymax a numeric value indicating the maximum y value plotted in each graph
#' @param xmin a numeric value indicating the minimum x value plotted in each graph
#' @param xmax a numeric value indicating the maximum x value plotted in each graph
#' @param xplots a numeric value indicating how many plots to plot in each row
#' @param height a numeric value indicating the height of the pdf file generated
#' @param width a numeric value indicating the width of the pdf file generated
#' @param cexPt a numeric value indicating the size to plot for points
#' @param cexX a numeric value indicating the size of x-axis text
#' @param cexY a numeric value indicating the size of y-axis text
#' @param nameVector either a logial value indicating whether to plot the photograph names above the graph or not or a vector the same length as the number of pictures containing the desired names. Defaults to TRUE.
#' @param plotStandardLoc a logical value indicating whether to draw a dashed horizontal line at the standardization point
#' @param showNum a logical value indicating whether to annotate each graph with the picture number (determined alphabetically from the picture names)
#' @param popUp a logical value indicating whether to pop up the figure after it has been created
#' @param overwrite a logical value indicating whether to overwrite existing figures created on the same day for the same project name
#' @param savePDF a logical value indicating whether to save a PDF file or open a new quartz. Defaults to TRUE.

#' @return A pdf file with one plot for each photograph is saved to visualize the results of imageJ analyses

#' @export

#' @author Aleeza c. Gerstein

plotRaw <- function(projectName, ymin = 0, ymax=250, xmin = 0, xmax = 40, xplots = 6, height =4, width = 8, cexPt = 0.6, cexX = 0.8, cexY = 0.8, nameVector = TRUE , plotStandardLoc=TRUE, showNum=FALSE, popUp = TRUE, overwrite=TRUE, savePDF= TRUE){
	dir.create(paste("figures/", projectName,  sep=""), showWarnings = FALSE)
	t <- file.path("figures", projectName,  paste(projectName, "_raw.pdf", sep=""))
	if (!overwrite){
		if (file.exists(t)){
			t <- file.path("figures", projectName, paste(projectName, "_raw_2.pdf", sep=""))
			if (file.exists(t)){
				k <- 2
				while(file.exists(t)){
					k <- k+1
					t <- file.path("figures", projectName, paste(projectName, "_raw_", k, ".pdf", sep=""))
					}
				}
			}
		}
	data <- eval(parse(text=projectName))
	standardLoc <- 2.5
	dotMax <- max(sapply(data, function(x) {x[which(x[,1] > standardLoc)[1], 2]})) 		
	standards <-c( sapply(data, function(x) {dotMax-x[which(x[,1] > standardLoc)[1], 2]}))	
	convert <- unlist(lapply(data, function(x) 40/length(x[,1])))
	if (is.logical(nameVector)){
		if (nameVector){label <- names(data)}		
		else {label <- rep("", length(data))}
		}
	else {label <- nameVector}

	if (xplots > length(data)){
		xplots <- length(data)
		}
	if (ceiling(length(data)/xplots) < 6) {
		yplots<- ceiling(length(data)/xplots)}
	else {yplots<- 6}
	numpages <- ceiling(length(data)/(xplots*yplots))
	if(savePDF){
		pdf(t, width=width, height=height)
		}
	# if(!savePDF){
		# quartz(width=width, height=height)
		# }
	par(mfrow=c(yplots , xplots), mar=c(1,1,1,1), oma=c(4,5,1,1))
	for (i in 1:length(data)){
		.discplotNoRep(data[[i]], label[i], ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax, stand=standards[i], standardLoc = standardLoc, cexPt = cexPt, plotStandardLoc = plotStandardLoc)
		if(numpages == 1){
			if (i >= xplots*yplots-xplots+1){
				axis(1, cex.axis=cexX, at=c(0, 10, 20, 30, 40), labels=c(0, 10, 20, 30, 40))
			}
		}
		if(numpages == 2){
			if (i >= xplots*yplots-xplots+1 & i < xplots*yplots+1){
				axis(1, cex.axis=cexX, at=c(0, 10, 20, 30, 40), labels=c(0, 10, 20, 30, 40))
			}
			if (i >= 2*xplots*yplots-xplots+1){
				axis(1, cex.axis=cexX, at=c(0, 10, 20, 30, 40), labels=c(0, 10, 20, 30, 40))
			}
		}				
		if(numpages == 3){
			if (i >= xplots*yplots-xplots+1 & i < xplots*yplots+1){
				axis(1, cex.axis=cexX, at=c(0, 10, 20, 30, 40), labels=c(0, 10, 20, 30, 40))
			}
			if (i >= 2*xplots*yplots-xplots+1 & i < 2*xplots*yplots+1){
				axis(1, cex.axis=cexX, at=c(0, 10, 20, 30, 40), labels=c(0, 10, 20, 30, 40))
			}
			if (i >= (length(data)-xplots)){
				axis(1, cex.axis=cexX, at=c(0, 10, 20, 30, 40), labels=c(0, 10, 20, 30, 40))
			}
		}				
		k <- 1
		while (k <= numpages){
		if (i %in% seq(1, k*yplots*xplots, by=xplots)) {axis(2, cex.axis=cexY, las=2)}
			k <- k+1}
		if(showNum){
			text(maxDist*0.95, maxDist*0.95, i)
		}
	}
	mtext("Distance from center of disk (mm)", side= 1, outer=TRUE, line=2)
	mtext("Pixel intensity", side=2, outer=TRUE, line=2)
	if(savePDF){
		dev.off()
		cat(paste("\tFigure saved: ", t, sep=""))
		if(popUp){
		tt <- paste("open ",t)
		system(tt)
	}
	}
	}

.discplotNoRep <- function(data,  label=label, ymin=0, ymax=250, xmin=0, xmax=40, standardLoc = 2.5, cexPt = 0.6, stand = 0, xaxt="n", yaxt="n", plotStandardLoc =FALSE){
	plot(data[,1], data[,2]+stand, ylim=c(ymin, ymax), xlim=c(xmin, xmax), xaxt=xaxt, yaxt="n", cex=cexPt, col="black")
	if (yaxt=="s"){
		axis(2, las=2)}
	axis(1, labels=FALSE, at=c(0, 10, 20, 30, 40))
	axis(2, labels=FALSE)
	mtext(label, side=3, cex=0.6)
	if(plotStandardLoc){
		abline(v= standardLoc, lty=2, col="red")
		}
	}
