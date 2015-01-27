#' Used to plot the ZOI, slope and AUC parameter results

#' @description This function creates a pdf figure of plots showing the results of the imageJ analysis for resistance (ZOI) and tolerance (AUC).

#' @inheritParams plotRaw
#' @inheritParams twoParamPlot

#' @details Basic parameter plotting functions for three parameter plots (ZOI, fAUC , slope). Input can be the dataframe from either \code{\link{createDataframe}} \code{type="df"} or from \code{\link{aggregateData}} \code{type=="ag"}. The default is to plot tolerance as a barplot and ZOI and slope as a dotplot, tolerance can also be plotted as a dotplot with \code{barplot=FALSE} though there is currently not support to plot either ZOI or slope as a barplot in this framework. 

#' @return Either a pdf figure figure saved to the 'figures' directory ("projectName_ZOI-slope-fAUC.pdf" or a figure on screen

#' @export

#' @author Aleeza C. Gerstein


threeParamPlot <- function(projectName, type, ZOI = "ZOI20", AUC = "fAUC20", ZOImin = 30, slopeMax = 160, tolMax = 100, width = 6, height = 4, xlabels="line", xlabAngle=NA, order=NA, orderFactor = "line", overwrite=TRUE, savePDF= TRUE, popUp = TRUE){
	
	dir.create(paste("figures/", projectName,  sep=""), showWarnings = FALSE)
	t <- file.path("figures", projectName,  paste(projectName, "_ZOI-slope-fAUC.pdf", sep=""))
	if (!overwrite){
		if (file.exists(t)){
			t <- file.path("figures", projectName, paste(projectName, "_ZOI-slope-fAUC_2.pdf", sep=""))
			if (file.exists(t)){
				k <- 2
				while(file.exists(t)){
					k <- k+1
					t <- file.path("figures", projectName, paste(projectName, "_ZOI-slope-fAUC_", k, ".pdf", sep=""))
					}
				}
			}
		}
	if(type == "ag" & !is.na(order[1])){
		data <- eval(parse(text=paste(projectName, ".ag", sep="")))	
		var <- substring(names(data)[length(data)], 1, 2)
		if(order[1]=="factor"){
			ordData<-data[order(data[, orderFactor]),] 
			if(length(xlabels)==1){
		 		xlabels <- as.character(ordData[, xlabels])
			}	
		}
		if(!order[1]=="factor"){
			ordData <-  data[order, ]
			if(length(xlabels)==1){
				 xlabels <- as.character(ordData[, xlabels])
			}
		}
		
	}
	if(is.na(order[1])){
		if(type=="ag"){
			ordData <- eval(parse(text=paste(projectName, ".ag", sep="")))	
			var <- substring(names(ordData)[length(ordData)], 1, 2)	
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
	par(mfrow=c(3, 1), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))
	if(type=="ag"){
		plot(mp[1,], ordData[, ZOI], ylim=c(ZOImin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
		arrows(mp[1,], ordData[, ZOI]-ordData[, paste(var, ".", ZOI, sep="")], mp[1,], ordData[, ZOI]+ordData[,paste(var, ".", ZOI, sep="")], length=0)
	axis(1, at=mp[1,], labels=FALSE)
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, ZOI], ylim=c(ZOImin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(xlabels)+0.5))
	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
	}
	axis(2, las=2, cex.axis=1)
	mtext("Distance from\n disk (mm)", side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(A), " Resistance", sep="")), side=3, adj=0.01)

	if(type=="ag"){
		plot(mp[1,], ordData[, "slope"], ylim=c(0, slopeMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
		arrows(mp[1,], ordData[, "slope"]-ordData[, paste(var, ".slope", sep="")], mp[1,], ordData[, ZOI]+ordData[,paste(var, ".slope", sep="")], length=0)
	axis(1, at=mp[1,], labels=FALSE)
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, "slope"], ylim=c(0, slopeMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(xlabels)+0.5))
	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
	}
	axis(2, las=2, cex.axis=1)
	title <- as.list(expression("slope" , paste( "(", Delta, "intensity/", Delta, "distance)", sep="")))
	mtext(do.call(expression, title), side=2, cex=0.8, line = c(3.75,2.5))
	mtext(expression(paste(bold(B), " Sensitivity", sep="")), side=3, adj=0.01)
		
	if(type=="ag"){	
		mp <- barplot(t(tols*100), ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
		box()
		 arrows(mp[1,], ordData[,AUC]*100-ordData[,paste(var, ".", AUC, sep="")]*100, mp[1,], ordData[,AUC]*100+ ordData[,paste(var, ".", AUC, sep="")]*100, length=0)
		if(is.na(xlabAngle)) 	axis(1, at=mp[1,], labels=xlabels)
		else{
			axis(1, at=mp[1,], labels=FALSE)
			text(mp[1,],  -10, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, AUC]*100, ylim=c(0, tolMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(unique(xlabels))+0.5))
		if(is.na(xlabAngle)){
			 axis(1, at=1:length(xlabels), labels=xlabels)
			 print("here")
			 }
		else{
			axis(1, at=1:length(xlabels), labels=FALSE)
			text(1:length(xlabels),  -5, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	axis(2, las=2, at=c(0, 20, 40, 60, 80, 100), cex.axis=1)
	mtext("Growth above\n ZOI (%)",  side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(C), " Tolerance", sep="")), side=3, adj=0.01)
	if(savePDF){
		dev.off()
		cat(paste("\tFigure saved: ", t, sep=""))
		if(popUp){
			tt <- paste("open ",t)
			system(tt)
		}
	}

}