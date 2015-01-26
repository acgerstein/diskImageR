#' Used to plot the ZOI and AUC results

#' @description This function creates a pdf figure of plots showing the results of the imageJ analysis for resistance (ZOI) and tolerance (AUC).

#' @inheritParams plotRaw
#' @param ZOI specify the ZOI parameter to be plotted ("ZOI20", "ZOI50" or "ZOI80"), default = "ZOI20".
#' @param AUC specify the AUC parameterto be plotted ("fAUC20", "fAUC50" or "fAUC80"), default = "fAUC20".
#' @param ZOImin minimum distance from the disk for resistance plot (minimum y axis value), default = 30.
#' @param ZOImin minimum distance from the disk for resistance plot (minimum y axis value), default = 30.
#' @param tolMax maximum y axis value for tolerance plot. Note tolerance is coverted to a perent, default = 100.
#' @param xlabels either a vector containing the desired x-axis labels, or a single value indicating the column name that contains the values to use (likely either the 'line' column or one of the type columns), default = "line".
#' @param xlabAngle indicates whether to print the x axis labels on a angle, if a number is provided this will be the angle used. The defauilt is not to plot on an angle, default = NA.
#' @param order can be either "factor" or "custom". If custom, supply a numberial vector the same length as the dataframe to indicate the desired order. If factor, supply the column name in \code{ordeFactor} to be used to factor. 
#' @param orderFactor if \code{order = "factor"} supply the column name to be used to factor. 
#' @param barplot whether to plot tolerance as a barplot (barplot = TRUE) or dotplot (barplot = FALSE), default = TRUE. Only possible when \code{type = "ag"}

#' @details Basic parameter plotting functions to plot ZOI and fAUC parameter plots. Input can be the dataframe from either \code{\link{createDataframe}} \code{type="df"} or from \code{\link{aggregateData}} \code{type=="ag"}. The default is to plot ZOI as a dotplot and tolerance as a barplot, though tolerance can also be plotted as a dotplot with \code{barplot=FALSE} (currently there is not support to plot ZOI as a barplot in this framework). 

#' @return Either a pdf figure figure (projectName_ZOI-fAUC.pdf) saved to the 'figures' directory or a figure on screen

#' @export

#' @author Aleeza C. Gerstein

twoParamPlot <- function(projectName, type, ZOI = "ZOI20", AUC = "fAUC20",  ZOImin = 30, tolMax = 100, width = 6, height = 4, xlabels ="line", xlabAngle=NA, order=NA, orderFactor = "line", overwrite=TRUE, savePDF= TRUE, popUp = TRUE, barplot=TRUE){
	if(!(hasArg(type))){
		cont <- readline(paste("Please select whether dataframe is from 'createDataframe' (df) or `aggregateData (ag) ", sep=""))
		type <- cont
	}
	print(orderFactor)
	print(xlabels)
	dir.create(paste("figures/", projectName,  sep=""), showWarnings = FALSE)
	t <- file.path("figures", projectName,  paste(projectName, "_ZOI-fAUC.pdf", sep=""))
	if (!overwrite){
		if (file.exists(t)){
			t <- file.path("figures", projectName, paste(projectName, "_ZOI-fAUC_2.pdf", sep=""))
			if (file.exists(t)){
				k <- 2
				while(file.exists(t)){
					k <- k+1
					t <- file.path("figures", projectName, paste(projectName, "_ZOI-fAUC_", k, ".pdf", sep=""))
					}
				}
			}
		}
	print(xlabels)
	if(type == "ag" & !is.na(order[1])){
		data <- eval(parse(text=paste(projectName, ".ag", sep="")))	
		var <- substring(names(data)[length(data)], 1, 2)
		print(xlabels)
		if(order[1]=="factor"){
			ordData<-data[order(data[, orderFactor]),] 
			if(length(xlabels)==1){
		 		xlabels <- as.character(ordData[, xlabels])
			}	
		}
		if(!order[1]=="factor"){
			ordData <-  data[order, ]
			if(length(xlabels)==1){
				print("here")
				 xlabels <- as.character(ordData[, xlabels])
				 xlabels <- xlabels[order, ]
				print(xlabels)
			}
		}
		
	}
	if(is.na(order[1])){
		if(type=="ag"){
			var <- substring(names(data)[length(data)], 1, 2)
			ordData <- eval(parse(text=paste(projectName, ".ag", sep="")))	
			if(length(xlabels)==1){
				 xlabels <- as.character(ordData[, xlabels])
			}
		}
		
		if(type=="df"){
			ordData <- eval(parse(text=paste(projectName, ".df", sep="")))
			if(length(xlabels)==1){
				 xlabels <- unique(as.character(ordData[, xlabels]))
			}
		}
	}
	tols <- ordData[, AUC]
	mp <- barplot(t(tols), beside=TRUE, plot=FALSE)	
	if(savePDF){
		 pdf(t, width=width, height=height)
		}	
	par(mfrow=c(2, 1), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))
	
	if(type=="ag"){
		if(barplot == TRUE){		
			plot(mp[1,], ordData[, ZOI], ylim=c(ZOImin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
			arrows(mp[1,], ordData[, ZOI]-ordData[, paste(var, ".", ZOI, sep="")], mp[1,], ordData[, ZOI]+ordData[,paste(var, ".", ZOI, sep="")], length=0)
		axis(1, at=mp[1,], labels=FALSE)
		}
	else{
			plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, ZOI], ylim=c(ZOImin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
			arrows(as.numeric(as.factor(ordData[, orderFactor])), ordData[, ZOI]-ordData[, paste(var, ".", ZOI, sep="")], as.numeric(as.factor(ordData[, orderFactor])), ordData[, ZOI]+ordData[,paste(var, ".", ZOI, sep="")], length=0)
		axis(1, at=mp[1,], labels=FALSE)
		}
		
	}
	
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, ZOI], ylim=c(ZOImin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(xlabels)+0.5))
	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
	}
	
	axis(2, las=2, cex.axis=0.8)
	mtext("Distance\n from disk (mm)", side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(A), " Resistance", sep="")), side=3, adj=0.01)
	
	if(type=="ag"){
		if(barplot == TRUE){	
			mp <- barplot(t(tols*100), ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
			box()
		 	arrows(mp[1,], ordData[,AUC]*100-ordData[, paste(var, ".", AUC, sep="")]*100, mp[1,], ordData[,AUC]*100+ ordData[,paste(var, ".", AUC, sep="")]*100, length=0)
			if(is.na(xlabAngle)) 	axis(1, at=mp[1,], labels=xlabels)
			else{
				axis(1, at=mp[1,], labels=FALSE)
				text(mp[1,],  -10, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
			}
		}
		else{
			plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, fAUC], ylim=c(0, tolMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(xlabels)+0.5))
			axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
			arrows(as.numeric(as.factor(ordData[, orderFactor])), ordData[,AUC]*100-ordData[, paste(var, ".", AUC, sep="")]*100, as.numeric(as.factor(ordData[, orderFactor])), ordData[,AUC]*100+ ordData[,paste(var, ".", AUC, sep="")]*100, length=0)
			if(is.na(xlabAngle)) 	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=xlabels)
			else{
				axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
				text(as.numeric(as.factor(unique(ordData[, orderFactor]))),  -10, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
			}
		 }
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, AUC]*100, ylim=c(0, tolMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(unique(xlabels))+0.5))
		if(is.na(xlabAngle)){
			 axis(1, at=1:length(xlabels), labels=xlabels)
			 }
		else{
			axis(1, at=1:length(xlabels), labels=FALSE)
			text(1:length(xlabels),  -5, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	axis(2, las=2, at=c(0, 20, 40, 60, 80, 100), cex.axis=0.8)
	mtext("Growth\n above ZOI (%)",  side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(B), " Tolerance", sep="")), side=3, adj=0.01)
	if(savePDF){
		dev.off()
		cat(paste("\tFigure saved: ", t, sep=""))
		if(popUp){
			tt <- paste("open ",t)
			system(tt)
		}
	}

}
