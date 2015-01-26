oneParamPlot <- function(projectName, type, param  = "ZOI20", ymin = 0, ymax = 100, width = 6, height = 4, xlabels="line", xlabAngle=NA, order=NA, orderFactor = "line", overwrite=TRUE, savePDF= TRUE, popUp = TRUE){
	
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

	if(type == "ag" & !is.na(order)){
		data <- eval(parse(text=paste(projectName, ".ag", sep="")))	
		var <- substring(names(data)[length(data)], nchar(names(data)[length(data)])-1, nchar(names(data)[length(data)]))
		if(order=="factor"){
			ordData <- arrange(data, data[, orderFactor])
			}
		if(order=="custom"){
			ordData <- data %>%
						 		mutate(order) %>%
									arrange(order)			
		}
	}
	else{
		if(type=="ag"){
			var <- substring(names(data)[length(data)], 1, 2)
			ordData <- eval(parse(text=paste(projectName, ".ag", sep="")))	
			}
		if(type=="df"){
			ordData <- eval(parse(text=paste(projectName, ".df", sep="")))
			}
	}
	tols <- ordData[, AUC]
	mp <- barplot(t(tols), beside=TRUE, plot=FALSE)	
	if(length(xlabels)==1){
		 print(xlabels)
		 xlabels <- unique(as.character(ordData[, xlabels]))
		}
	print(xlabels)
	if(grep("ZOI", param){
		yrange <- c(ymax, ymin)
		}
	else{
		yrange <- c(ymin, ymax)
		}
	if(savePDF){
		 pdf(t, width=width, height=height)
		}
	par(mfrow=c(2, 1), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))

	if(type=="ag"){
			plot(mp[1,], ordData[, param], ylim=yrange, yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
				arrows(mp[1,], ordData[, ZOI]-ordData[, paste(var, ".", ZOI, sep="")], mp[1,], ordData[, ZOI]+ordData[,paste(var, ".", ZOI, sep="")], length=0)
		axis(1, at=mp[1,], labels=FALSE)
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, ZOI], ylim=yrange, yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(unique(ordData[, orderFactor]))+0.5))
	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
	}
	axis(2, las=2, cex.axis=0.8)
	mtext("Distance\n from disk (mm)", side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(A), " Resistance", sep="")), side=3, adj=0.01)
	
	if(type=="ag"){	
		mp <- barplot(t(tols*100), ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
		box()
		 arrows(mp[1,], ordData[,AUC]*100-ordData[,paste(var, ".", AUC, sep="")]*100, mp[1,], ordData[,AUC]*100+ ordData[,paste(var, ".", AUC, sep="")]*100, length=0)
		if(is.na(xlabAngle)) 	axis(1, at=mp[1,], labels=xlabels)
		else{
			axis(1, at=mp[1,], labels=FALSE)
			text(mp[1,],  -5, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, AUC]*100, ylim=c(0, tolMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(unique(ordData[, orderFactor]))+0.5))
		if(is.na(xlabAngle)){
			 axis(1, at=1:length(xlabels), labels=xlabels)
			 print("here")
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

