#' Used to plot a single parameter

#' @description This function creates a pdf figure of plots showing the results of the imageJ analysis for resistance (ZOI), sensitivity (slope) and tolerance (AUC).

#' @inheritParams plotRaw
#' @inheritParams twoParamPlot
#' @param barplot whether to plot values as a barplot (barplot = TRUE) or dotplot (barplot = FALSE), default = TRUE. Only possible when \code{type = "ag"}

#' @details Basic parameter plotting functions to plot a single  parameter. Input can be the dataframe from either \code{\link{createDataframe}} \code{type="df"} or from \code{\lnk{aggregateData}} \code{type=="ag"}. 

#' @return Either a pdf figure figure (projectName_ZOI-fAUC.pdf) saved to the 'figures' directory or a figure on screen

#' @export

#' @author Aleeza C. Gerstein

oneParamPlot <- function(projectName, type, param  = "ZOI20", ymin = 0, ymax = 100, width = 6, height = 4, xlabels="line", xlabAngle=NA, order=NA, orderFactor = "line", overwrite=TRUE, savePDF= TRUE, popUp = TRUE, barplot = TRUE){
	
	dir.create(paste("figures/", projectName,  sep=""), showWarnings = FALSE)
	t <- file.path("figures", projectName,  paste(projectName, "_", param, ".pdf", sep=""))
	if (!overwrite){
		if (file.exists(t)){
			t <- file.path("figures", projectName, paste(projectName, "_", param, "_2.pdf", sep=""))
			if (file.exists(t)){
				k <- 2
				while(file.exists(t)){
					k <- k+1
					t <- file.path("figures", projectName, paste(projectName, "_", param, "_", k, ".pdf", sep=""))
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



	if(length(xlabels)==1){
		 print(xlabels)
		 xlabels <- unique(as.character(ordData[, xlabels]))
		}
	print(xlabels)
	if(grep("ZOI", param)){
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
		if(barplot!=TRUE){
			plot(mp[1,], ordData[, param], ylim=yrange, yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
			arrows(mp[1,], ordData[, ZOI]-ordData[, paste(var, ".", ZOI, sep="")], mp[1,], ordData[, ZOI]+ordData[,paste(var, ".", ZOI, sep="")], length=0)
		axis(1, at=mp[1,], labels=FALSE)
		}
	else{
		tols <- ordData[, param]
		if(!param %in% c("fAUC20", "fAUC50", "fAUC80")){
			mp <- barplot(t(tols), ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
			}
		else{
			mp <- barplot(t(tols)*100, ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
			}
		arrows(mp[1,], ordData[,AUC]*100-ordData[,paste(var, ".", AUC, sep="")]*100, mp[1,], ordData[,AUC]*100+ ordData[,paste(var, ".", AUC, sep="")]*100, length=0)
		if(is.na(xlabAngle)) 	axis(1, at=mp[1,], labels=xlabels)
		else{
			axis(1, at=mp[1,], labels=FALSE)
			text(mp[1,],  -5, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, ZOI], ylim=yrange, yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(unique(ordData[, orderFactor]))+0.5))
	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
	}
	mtext("Distance\n from disk (mm)", side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(A), " Resistance", sep="")), side=3, adj=0.01)
	if(is.na(xlabAngle)){
		 axis(1, at=1:length(xlabels), labels=xlabels)
		 print("here")
		 }
	else{
		axis(1, at=1:length(xlabels), labels=FALSE)
		text(1:length(xlabels),  -5, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	axis(2, las=2, cex.axis=0.8)
	if(param %in% c("fAUC20", "fAUC50", "fAUC80")) mtext("Growth\n above ZOI (%)",  side=2, line=2.5, cex=0.8)
	if(param %in% c("ZOI20", "ZOI50", "ZOI80")) mtext("Distance from\n disk (mm)", side=2, line=2.5, cex=0.8)
	if(param=="slope"){
		title <- as.list(expression("slope" , paste( "(", Delta, "intensity/", Delta, "distance)", sep="")))
		mtext(do.call(expression, title), side=2, cex=0.8, line = c(3.75,2.5))
		}

	if(savePDF){
		dev.off()
		cat(paste("\tFigure saved: ", t, sep=""))
		if(popUp){
			tt <- paste("open ",t)
			system(tt)
		}
	}

}
