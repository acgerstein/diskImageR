#' Used to plot the parameter results 

#' @description This function creates a pdf figure of plots showing the results of the imageJ analysis, with one plot for each picture. This function is primarily for visualization purposes.

#' @inheritParams plotRaw
#' @inheritParams plotRaw
#' @param ZOImin
#' @param tolMax
#' @param slopeMax
#' @param xlabels
#' @param xlabAngle
#' @param order

#' @deails

#' @return Either a pdf figure figure saved to the 'figures' directory or a figure on screen

#' @export

#' @author Aleeza c. Gerstein

twoParamPlot <- function(projectName, ZOI = "ZOI20", AUC = "fAUC20", var = "se", ZOImin = 30, tolMax = 100, slopeMax = 200, width = 6, height = 4, xlabels=mp[1,], xlabAngle=NA, order=NA, orderFactor = NA overwrite=TRUE){
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

	data <- eval(parse(text=projectName))
	type <- substring(t, nchar(t)-1, nchar(t))
	if(type == "ag"){
		if(!is.na(order)){
			if(order=="factor"){
				ordData <- arrange(data, data[, orderFactor])
				}
			if(order=="custom"){
				ordData <- data %>%
							 		mutate(order) %>%
										arrange(order)			
			}
		}
		else {
			tols <- data[, AUC]
			mutate(data, 1:length(data[, AUC]))
			}
		tols <- ordData[, AUC]
		mp <- barplot(t(tols), beside=TRUE, plot=FALSE)	
	}
	
	pdf(t, width=width, height=height)
	par(mfrow=c(2, 1), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))
	if(type=="ag"){
		plot(mp[1,], ordData[, ZOI], ylim=c(ZOImin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
		arrows(mp[1,], ordData[, ZOI]-ordData[,paste(var, ".", ZOI, sep=""), mp[1,], ordData[, ZOI]+ordData[,paste(var, ".", ZOI, sep=""), length=0)
	axis(1, at=mp[1,], labels=FALSE)
	}
	if(type=="df"){
		plot(ordData[, orderFactor], ordData[, ZOI], ylim=c(ZOImin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4)
	axis(1, at=as.numeric(unique(ordData[, orderFactor]), labels=FALSE)
	}
	axis(2, las=2, cex.axis=0.8)
	mtext("Distance from disk (mm)", side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(A), " Resistance", sep="")), side=3, adj=0.01)
	
	if(type=="ag"){	
		mp <- barplot(t(tols*100), ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
		box()
		arrows(mp[1,], ordData[,AUC]*100-ordData[,paste(var, ".", AUC, sep="")*100, mp[1,], ordData[,AUC]*100+ ordData[,paste(var, ".", AUC, sep="")*100, length=0)
		if(is.na(labAngle)) 	axis(1, at=mp[1,], labels=xlabels)
		else{
			axis(1, at=mp[1,], labels=FALSE)
			text(mp[1,],  -5, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	if(type=="df"){
		plot(as.numeric(ordData[, orderFactor]), ordData[, AUC], ylim=c(0, tolMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4)
		if(is.na(labAngle)) 	axis(1, at=as.numeric(unique(ordData[, orderFactor])), labels=xlabels)
		else{
			axis(1, at=as.numeric(unique(ordData[, orderFactor])), labels=FALSE)
			text(unique(ordData[, orderFactor])),  -5, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	axis(2, las=2, at=c(0, 20, 40, 60, 80, 100), cex.axis=0.8)
	mtext("Growth above ZOI (%)",  side=2, line=2.5, cex=0.8)
		mtext(expression(paste(bold(B), " Tolerance", sep="")), side=3, adj=0.01)
	dev.off()	
	system(paste("open ", figureName, sep=""))
}

# threeParamPlot <- function(dat.ag, order=NA, figureName, ZOI = "ZOI20", AUC = "fAUC20", var = "se", ZOImin = 60, slopeMax = 200, tolMax = 100, xlabels=mp[1,], plotLegend=FALSE, legendText="", addLines=FALSE, labAngle=NA, width = 5, height = 5){
	# if(!is.na(order)){
		# tols <- dat.ag[, AUC][order]	
		# }
	# else {
		# tols <- dat.ag[, AUC]
		# order <- 1:length(dat.ag[,1])
		# }
	# mp <- barplot(t(tols), beside=TRUE, plot=FALSE)	
	# pdf(figureName, width=width, height=height)
	# par(mfrow=c(3, 1), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))
	# plot(mp[1,], dat.ag[, ZOI][order], ylim=c(ZOImin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), xlim=c(0, max(mp)+1), cex=1.4)
	# axis(2, las=2, cex.axis=0.8)
	# arrows(mp[1,], dat.ag[, ZOI][order]-dat.ag[,paste(var, ".", ZOI, sep="")][order], mp[1,], dat.ag[, ZOI][order]+dat.ag[,paste(var, ".", ZOI, sep="")][order], length=0)
	# if(plotLegend){
		# legend(legendText)
		# }
	# if(addLines){
		# abline(h=14, lty=2)
		# abline(h=19, lty=2)
		# }
	# mtext("Distance (mm)", side=2, line=2.5, cex=0.8)
	# mtext(expression(paste(bold(A), " Resistance", sep="")), side=3, adj=0.01)
	# axis(1, at=mp[1,], labels=FALSE)

	# plot(mp[1,], dat.ag[, "slope"][order], ylim=c(0, slopeMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), xlim=c(0, max(mp)+1), cex=1.4)
	# axis(2, las=2, cex.axis=0.8)
	# arrows(mp[1,], dat.ag[, "slope"][order]-dat.ag[,paste(var, ".slope", sep="")][order], mp[1,], dat.ag[, "slope"][order]+dat.ag[,paste(var, ".slope", sep="")][order], length=0)
	# if(plotLegend){
		# legend(legendText)
		# }
	# if(addLines){
		# abline(h=14, lty=2)
		# abline(h=19, lty=2)
		# }
	# mtext(expression(paste("(", Delta, "intensity/", Delta, "distance)", sep="")), side=2, line=2.5, cex=0.8)
	# mtext("slope", side=2, line=4, cex=0.8)
	# mtext(expression(paste(bold(B), " Sensitivity", sep="")), side=3, adj=0.01)
	# axis(1, at=mp[1,], labels=FALSE)
	
	# mp <- barplot(t(tols*100), ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
	# box()
	# if(is.na(labAngle)) 	axis(1, at=mp[1,], labels=xlabels)
	# else{
			# axis(1, at=mp[1,], labels=FALSE)
			# text(mp[1,],  -5, xlabels, srt=labAngle, xpd=NA, adj=0, cex=0.8)
			# }
	# axis(2, las=2, at=c(0, 20, 40, 60, 80, 100), cex.axis=0.8)
	# mtext("Growth above ZOI (%)",  side=2, line=2.5, cex=0.8)
	# arrows(mp[1,], dat.ag[,AUC][order]*100-dat.ag[,paste(var, ".", AUC, sep="")][order]*100, mp[1,], dat.ag[,AUC][order]*100+dat.ag[,paste(var, ".", AUC, sep="")][order]*100, length=0)
	# mtext(expression(paste(bold(C), " Tolerance", sep="")), side=3, adj=0.01)
	# dev.off()	
	# system(paste("open ", figureName, sep=""))
# }

# basic2panel2means <- function(dat.ag, dat.ag2, order, figureName, ZOI = "ZOI10", AUC = "fracAUC10", ZOImin = 60, tolMax = 100, xlabels=mp[1,], plotLegend=FALSE, legendText="", addLines=FALSE){
	# tols <- cbind(dat.ag[, AUC][order], dat.ag2[,AUC][order])
	# mp <- barplot(t(tols), beside=TRUE, plot=FALSE)	
	# pdf(figureName, width=6.5, height=6)
	# par(mfrow=c(2, 1), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))
	# plot(mp[1,], dat.ag[, ZOI][order], ylim=c(ZOImin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), xlim=c(0, max(mp)+1), cex=1.4)
	# points(mp[2,], dat.ag2[,ZOI][order], col=grey(0.7), pch=19, cex=1.4)
	# axis(2, las=2, cex.axis=0.8)
	# arrows(mp[1,], dat.ag[, ZOI][order]-dat.ag[,paste("se.", ZOI, sep="")][order], mp[1,], dat.ag[, ZOI][order]+dat.ag[,paste("se.", ZOI, sep="")][order], length=0)
	# arrows(mp[2,], dat.ag2[, ZOI][order]-dat.ag2[,paste("se.", ZOI, sep="")][order], mp[2,], dat.ag2[, ZOI][order]+dat.ag2[,paste("se.", ZOI, sep="")][order], length=0)
	# if(plotLegend){
		# legend(legendText)
		# }
	# if(addLines){
		# abline(h=14, lty=2)
		# abline(h=19, lty=2)
		# }
	# mtext("Distance from disk (mm)", side=2, line=2.5)
	# mtext(expression(paste(bold(A), " Resistance", sep="")), side=3, adj=0.01)
	# axis(1, at=xlabels, labels=FALSE)
	
	# mp <- barplot(t(tols*100), ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
	# box()
	# axis(1, at=xlabels, labels=FALSE)
	# text(xlabels, -5, paste("A", order, sep=""), srt=-45, xpd=NA, adj=0, cex=0.8)
	# axis(2, las=2, at=c(0, 20, 40, 60, 80, 100), cex.axis=0.8)
	# mtext("Growth within ZOI (%)",  side=2, line=2.5)
	# arrows(mp[1,], dat.ag[,AUC][order]*100-dat.ag[,paste("se.", AUC, sep="")][order]*100, mp[1,], dat.ag[,AUC][order]*100+dat.ag[,paste("se.", AUC, sep="")][order]*100, length=0)
	# arrows(mp[2,], dat.ag2[,AUC][order]*100-dat.ag2[,paste("se.", AUC, sep="")][order]*100, mp[2,], dat.ag2[,AUC][order]*100+dat.ag2[,paste("se.", AUC, sep="")][order]*100, length=0)
	# mtext(expression(paste(bold(B), " Tolerance", sep="")), side=3, adj=0.01)
	# dev.off()	
	# system(paste("open ", figureName, sep=""))
# }

