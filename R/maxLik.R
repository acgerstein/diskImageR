#' Maximum Likelihood Inference

#' @description Find the maximum likelihood point of a single and double logistic equation by nonlinear optimisation to parameterize the single and double logistic models to use for determination of resistance, tolerance, and sensitivity

#' @param projectName the short name in use for the project.
#' @param clearHalo numeric value that indicates which picture should be used to represent a clear halo (i.e., the clear space beside the disk). This is extremely important to determine tolerance, as the intensity beside the disk for the photograph chosen is subtracted for all photographs. Choosing the correct photograph to be use is the only subjective aspect of this pipeline; lighting and camera settings will determine the degre to which the hue of the plate backbground changes among different photographs. Care should be taken to ensure that the plate background will be as similar as possible between plates, and thus any plate with a clear area beside the disk can be chosen. Photos are numbered alphabetically by name, and can also be found using \code{\link{plotRaw}}, showNum = TRUE. 
#' @param diskDiam the diameter of the diffusion disk in mm, defaults to 6.
#' @param maxDist a numeric value indicating the maximum distance away from the disk to be considered. Defaults to 30mm.
#' @param xplots a numeric value indicating how many plots to plot in each row, does not influence maximum likelihood fitting
#' @param ymax a numeric value indicating the maximum y value plotted in each graph, does not influence maximum likelihood fitting
#' @param height a numeric value indicating the height of the pdf file generated, does not influence maximum likelihood fitting
#' @param width a numeric value indicating the width of the pdf file generated, does not influence maximum likelihood fitting
#' @param AUC a numeric value indicating the the critical level of the area under the curve to plot, does not influence maximum likelihood fitting. Current only \code{AUC} = "80" (80\% reduction in growth), \code{AUC} = "50" (50\% reduction in growth), and \code{AUC} = "20" (20\% reduction in growth) are supported.
#' @param ZOI a numeric value indicating the the critical level of the zone of inhibition (resistance) parameter to plot, does not influence maximum likelihood fitting. Currently only \code{ZOI} = "80" (80\% reduction in growth), \code{ZOI} = "50" (50\% reduction in growth), \code{ZOI} = "20" (20\% reduction in growth), and \code{ZOI} = "all" are supported.
#' @param needML a logical value indicating whether the maximum likelihood results already exist in the global environment or not. If \code{\link{maxLik}} has already been run in this session then needML can be set to FALSE, which allows the user to replot the results without the need to rerun the time consuming maximum likelihood models. Defaults to TRUE.
#' @param popUp a logical value indicating whether to pop up the figure after it has been created.
#' @param nameVector either a logial value indicating whether to plot the photograph names above the graph or not or a vector the same length as the number of pictures containing the desired names. Defaults to TRUE.
#' @param overwrite a logical value indicating whether to overwrite existing figures created on the same day for the same project name.defaults to TRUE.
#' @param plotAUC a logical value indicating whether to plot the AUC or not. Defaults to TRUE
#' @param savePDF a logical value indicating whether to save a PDF file or open a new quartz. Defaults to TRUE.
#' @param plotSub allows you to plot only a subset of photographs - indicate with a vector the corresponding numeric indices of the data you wish to plot. Photographs are numbered alphabetically by name, and the photograph numbers can also be found by using the showNum option in \code{\link{plotRaw}}. Defaults to NA, which will plot data from all photographs. Note this does not affect the analysis component, all data is always analyzed.

#' @details \code{\link{maxLik}} searches for the maximum likelihood (ML) parameter for a single logistic and double logistic equation using the pixel intensity information previously determined from \code{\link{IJMacro}}. The equations fit are
#'  \deqn{single ('ML'): y = asymA*exp(scalA(x-od50A))\(1+exp(scalA(x-od50A)))+N(0, \sigma)}
#'  \deqn{double ('ML2'): y = asymA*exp(scalA(x-od50A))\(1+exp(scalA(x-od50A)))+asymB*exp(scalB(x-od50B)))\(1+exp(scalB(x-od50B)))+N(0, \sigma)}
#' where asymA and asymB are the two asymptotes, od50A and odB are the midpoints (of the two curves), scalA and scalB are the slopes at odA and odB divided by asymA/4 and asymB/4, respectively. Specifically, \code{\link{maxLik}} uses the\code{\link[subplex]{subplex}} method of \code{\link{optim}}, as implemented in \code{\link[diversitree]{find.mle}}.  The single logistic is the essentially the same as the double, yet fits only a single asymptote, midpoint and slope. The results from the double logistic fit are used in \code{\link{createDataframe}} to find the resistance points as well as to fit the area under the curve and thus tolerance, the single logistic is used to determine the midpoint of the curve which is later used to find sensitivity, i.e., the slope at this midpoint.

#' @return Two lists, ML and ML2 are saved to the global environment. A pdf file with one plot for each photograph is saved to visualize the results of curve fitting, zone of inhibition (resistance) and the area under the curve (tolerance).

#' @export

#' @author Aleeza C. Gerstein

# See if can do clearHalo automatically somehow

maxLik <- function(projectName, clearHalo, diskDiam = 6, maxDist=30, ymax=125, xplots = 5, height = 8,  width = 8, AUC=20, ZOI="all", needML = TRUE, popUp = TRUE, nameVector=TRUE, overwrite = TRUE, plotAUC = TRUE, savePDF= TRUE, plotSub = NA){
	options(warn=-1)
	if(!(hasArg(clearHalo))){
		cont <- readline(paste("Please specify photograph number with a clear halo: ", sep=""))
		clearHalo <- as.numeric(cont)
	}
	if(!AUC %in% c(80, 50, 20)){
		stop("Current suppported AUC values = 80, 50, 20")
		}
	if(!ZOI %in% c(80, 50, 20, "all")){
		stop("Current suppported ZOI values = 80, 50, 20")
		}
	fileFolder <- projectName
	dir.create(file.path(getwd(), "figures"), showWarnings= FALSE)
	dir.create(file.path(getwd(), "figures", fileFolder), showWarnings= FALSE)
	
	data <- eval(parse(text=projectName))
	standardLoc <- 2.5
	if (is.logical(nameVector)){
		if (nameVector){label <- names(data)}		
		else {label <- rep("", length(data))}
		}
	else {label <- nameVector}	

	if (!is.logical(standardLoc)){
		dotMax <- max(sapply(data, function(x) {x[which(x[,1] > standardLoc)[1], 2]})) 		
		stand <-c( sapply(data, function(x) {dotMax-x[which(x[,1] > standardLoc)[1], 2]}))
		}
	else{
		stand <- rep(0, length(data))
		}
	dotedge <- diskDiam/2+0.4
	if(needML){		
		ML <-lapply(c(1:length(data)), .getstatsLog, data=data, dotedge=dotedge, maxDist=maxDist, stand=stand, maxSlope=20)
		assign(paste(projectName, ".ML", sep=""), ML, envir=globalenv())
		cat(paste("\n", projectName, ".ML has been written to the global environment\n", sep=""))
		ML2 <- lapply(c(1:length(data)), .getstats2Log, data=data, dotedge=dotedge, maxDist=maxDist, stand=stand, maxSlope=20)
		assign(paste(projectName, ".ML2", sep=""), ML2, envir=globalenv())
		cat(paste("\n", projectName, ".ML2 has been written to the global environment\n", sep=""))
	}
	if(!needML){		
		MLt <- paste(projectName, ".ML", sep="") 
		MLt2 <- paste(projectName, ".ML2", sep="") 
		ML <- eval(parse(text=MLt))
		ML2 <- eval(parse(text=MLt2))
		cat(paste("\nUsing existing ML results ", MLt, " & ", MLt2, sep=""))		
		}

	if(plotAUC){
		clearHaloData <- data[[clearHalo]]
		startX <- which(clearHaloData[,1] > dotedge+0.5)[1]
		stopX <- which(clearHaloData[,1] > maxDist - 0.5)[1]
		clearHaloData <- clearHaloData[startX:stopX, 1:2]
		clearHaloData$x <- clearHaloData$x + stand[clearHalo] 
		clearHaloData$distance <- clearHaloData$distance - (dotedge+0.5)
		clearHaloStand <- clearHaloData[1,2]
		 			
		.plotAUC(projectName, ML=ML, ML2=ML2, dotedge = dotedge, stand = stand, standardLoc = standardLoc, maxDist = maxDist, ymax = ymax, clearHaloStand = clearHaloStand, AUC=AUC, ZOI=ZOI, height = height, width=width, xplots = xplots,label=label, overwrite = overwrite, popUp = popUp, plotAUC = plotAUC, savePDF = savePDF, plotSub = plotSub)
	}
	alarm()
}

.curve <-  function(asym, ic50,scal, x) {asym*exp(scal*(x-ic50))/(1+exp(scal*(x-ic50)))}

.curve2 <- function(asym, od50, scal, asymB, od50B, scalB, x) { asym*exp(scal*(x-od50))/(1+exp(scal*(x-od50)))+asymB*exp(scalB*(x-od50B))/(1+exp(scalB*(x-od50B)))} 

.getstatsLog <- function(i, data, stand, dotedge=dotedge, maxDist=maxDist, maxSlope=100){
	cat(".")
	startX <- which(data[[i]][,1] > dotedge+0.5)[1]
	stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
	data[[i]] <- data[[i]][startX:stopX, 1:2]
	data[[i]] <- subset(data[[i]], x != "NA")
	data[[i]]$x <- data[[i]]$x+ stand[i] -min(data[[i]]$x+stand[i])  #the micel only fits when it goes down to 0
	data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
	data[[i]]$distance <- log(data[[i]]$distance)
	sumsquares.fit <- function(theta){
		asym<-theta[[1]]
		ic50<-theta[[2]]
		scal<-theta[[3]]
		sigma<-theta[[4]]
		y<-data[[i]]$x
		x<-data[[i]]$distance
		res <- dnorm(y, (asym*exp(scal*(x-ic50))/(1+exp(scal*(x-ic50)))), sigma, log= T)
		sum(res)
	}
	lowOD <- min(data[[i]]$x)
	highOD <- quantile(data[[i]]$x, 0.99)
	lower <- c(highOD*0.8, 0, 0,0)
	upper <- c(highOD, max(data[[i]]$distance), maxSlope,maxSlope)
		
	par.tryA <-c(asym = 0.9*highOD, ic50 = log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2)
	par.tryB<-c(asym = 0.9*highOD, ic50 = log(maxDist)/4, scal = maxSlope*0.1, sigma = 0.2)
	par.tryC<-c(asym = 0.9*highOD, ic50 = log(maxDist)/2, scal =  maxSlope*0.01, sigma = 0.1)
	par.tryD<-c(asym = 0.9*highOD, ic50 = log(maxDist)/2, scal = maxSlope*0.1, sigma = 0.1)		

	mlpoint<-c()
	mlpointA<-diversitree::find.mle(sumsquares.fit,par.tryA, method="subplex",upper=upper,lower=lower,control=list(maxit=50000))
	mlpointB<-diversitree::find.mle(sumsquares.fit,par.tryB,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
	mlpointC<-diversitree::find.mle(sumsquares.fit,par.tryC,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
	mlpointD<-diversitree::find.mle(sumsquares.fit,par.tryD,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))

	mlpoint <- if (mlpointA$lnLik>mlpointB$lnLik) mlpointA else mlpointB
	mlpoint <- if (mlpointC$lnLik>mlpoint$lnLik) mlpointC else mlpoint
	mlpoint <- if (mlpointD$lnLik>mlpoint$lnLik) mlpointD else mlpoint
	mlpoint
}

.getstats2Log <- function(i, data, stand, dotedge=dotedge, maxDist=maxDist, maxSlope=100){
	cat(".")
	startX <- which(data[[i]][,1] > dotedge+0.5)[1]
	stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
	data[[i]] <- data[[i]][startX:stopX, 1:2]
	data[[i]] <- subset(data[[i]], x != "NA")
	data[[i]]$x <- data[[i]]$x+ stand[i] -min(data[[i]]$x+stand[i])  #the micel only fits when it goes down to 0
	data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
	data[[i]]$distance <- log(data[[i]]$distance)
	sumsquares.fit <- function(theta){
		asym<-theta[[1]]
		od50<-theta[[2]]
		scal<-theta[[3]]
		sigma<-theta[[4]]
		asymB<-theta[[5]]
		od50B<-theta[[6]]
		scalB<-theta[[7]]
		y<-data[[i]]$x
		x<-data[[i]]$distance
		res <- dnorm(y, (asym*exp(scal*(x-od50))/(1+exp(scal*(x-od50)))+asymB*exp(scalB*(x-od50B))/(1+exp(scalB*(x-od50B)))), sigma, log= T)
		sum(res)
	}
	lowOD <- min(data[[i]]$x)
	highOD <- quantile(data[[i]]$x, 0.99)
	lower <- c(0, 0, 0,0, 0, 0, 0)
	upper <- c(highOD, log(maxDist), maxSlope, 10, highOD,  log(maxDist), maxSlope)
	
	#This is conservative, keeping them symmetric 	
	par.tryA <-c(asym = 0.9*highOD, od50 = log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2, asymB = 0.9*highOD, od50B = log(maxDist)/4, scalB = maxSlope*0.01)
	par.tryB <-c(asym = 0.9*highOD, od50 = log(maxDist)/4, scal = maxSlope*0.1, sigma =  0.2, asymB = 0.9*highOD, od50B = log(maxDist)/4, scalB = maxSlope*0.1)
	par.tryC<-c(asym = 0.9*highOD, od50 = log(maxDist)/2, scal =  maxSlope*0.01, sigma = 0.1, asymB = 0.9*highOD,od50B = log(maxDist)/2, scal =  maxSlope*0.01)
	par.tryD<-c(asym = 0.9*highOD, od50 = log(maxDist)/2, scal =  maxSlope*0.1, sigma = 0.1, asymB = 0.9*highOD,od50B = log(maxDist)/2, scalB =  maxSlope*0.1)
	#Change asym and od50 
	par.tryE <-c(asym = 0.5*highOD, od50 =  log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2, asymB = 0.7*highOD, od50B =  log(maxDist)/2, scalB = maxSlope*0.01)
	par.tryF <-c(asym = 0.5*highOD, od50 =  log(maxDist)/4, scal = maxSlope*0.1, sigma =  0.2, asymB = 0.7*highOD, od50B =  log(maxDist)/2, scalB = maxSlope*0.01)
	par.tryG <-c(asym = 0.5*highOD, od50 =  log(maxDist)/4, scal = maxSlope*0.1, sigma =  0.2, asymB = 0.7*highOD, od50B =  log(maxDist)/2, scalB = maxSlope*0.1)
	par.tryH <-c(asym = 0.5*highOD, od50 =  log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2, asymB = 0.7*highOD, od50B =  log(maxDist)/2, scalB = maxSlope*0.01)
	
	mlpoint<-c()
	mlpointA<-diversitree::find.mle(sumsquares.fit,par.tryA, method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
	mlpointB<-diversitree::find.mle(sumsquares.fit,par.tryB,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
	mlpointC<-diversitree::find.mle(sumsquares.fit,par.tryC,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
	mlpointD<-diversitree::find.mle(sumsquares.fit,par.tryD,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
	mlpointE<-diversitree::find.mle(sumsquares.fit,par.tryE,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
	mlpointF<-diversitree::find.mle(sumsquares.fit,par.tryF,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
	mlpointG<-diversitree::find.mle(sumsquares.fit,par.tryG,method="subplex",upper=upper,lower=lower,control=list(maxit=50000))	
	mlpointH<-diversitree::find.mle(sumsquares.fit,par.tryH,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))	

	mlpoint <- if (mlpointA$lnLik>mlpointB$lnLik) mlpointA else mlpointB
	mlpoint <- if (mlpointC$lnLik>mlpoint$lnLik) mlpointC else mlpoint
	mlpoint <- if (mlpointD$lnLik>mlpoint$lnLik) mlpointD else mlpoint
	mlpoint <- if (mlpointE$lnLik>mlpoint$lnLik) mlpointE else mlpoint
	mlpoint <- if (mlpointF$lnLik>mlpoint$lnLik) mlpointF else mlpoint
	mlpoint <- if (mlpointG$lnLik>mlpoint$lnLik) mlpointG else mlpoint			
	mlpoint <- if (mlpointH$lnLik>mlpoint$lnLik) mlpointH else mlpoint			
	mlpoint
}

.singleAUC <- function(data, ML, ML2, stand, clearHaloStand, dotedge = 3.4, maxDist = 40, ymax = 200, AUC=50, ZOI=50, i, label, plotAUC = TRUE, showIC = TRUE){
	startX <- which(data[[i]][,1] > dotedge+0.5)[1]
	stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
	data[[i]] <- data[[i]][startX:stopX, 1:2]
	data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand 
	data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
	xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200) 	
	yy2.1<- .curve(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3],xx)
	yy2.2<- .curve(ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7],xx)		
	yy<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx) 
	#ZOI
	ploty <- data[[i]]$x
	ploty[ploty < 0] <-0
	slope <- ML[[i]]$par[3]
	ic50 <- ML[[i]]$par[2]	
	asym <- (ML[[i]]$par[1]+min(data[[i]]$x))
	plot(data[[i]]$distance, ploty, cex=0.7, col=grey(0.7), type="p", ylim=c(0, ymax), xlim=c(0, maxDist -dotedge), xaxt="n", yaxt="n", xlab="", ylab="")	
	axis(2, labels=FALSE)
	yyplot <- (yy+min(data[[i]]$x))
	yyplot[yyplot < 0] <- 0
	points(exp(xx), yyplot, type="l", col="black", lwd=3)			

	useAsym <- "TRUE"
	yy80halo <- yyplot[which.max(yyplot> asym * 0.2)]
	yy50halo <- yyplot[which.max(yyplot> asym * 0.5)]
	yy20halo <- yyplot[which.max(yyplot> asym * 0.8)]
	if(yy20halo < yy50halo){
		 yy20halo <- yyplot[which.max(yyplot> yyplot[length(yyplot)] * 0.8)]
		useAsym <- "FALSE"
	}

	xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200) 	
	xx80 <- exp(xx[which.max(yyplot> asym * 0.2)])
	xx50 <- exp(xx[which.max(yyplot> asym * 0.5)])
	xx20 <- exp(xx[which.max(yyplot> asym * 0.8)])
	if(useAsym == "FALSE"){
		 xx20 <- exp(xx[which.max(yyplot> yyplot[length(yyplot)] * 0.8)])
	}

	if(AUC==50){
		xx <- exp(xx[1:which.max(exp(xx) > xx50)-1])	
		}
	if(AUC==80){
		xx <- exp(xx[1:which.max(exp(xx) > xx80)-1])	
		}
	if(AUC==20){
		xx <- exp(xx[1:which.max(exp(xx) > xx20)-1])	
		}

	if(length(xx)<1){
		xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
	}		
	
	yy<- .curve2(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], log(xx)) 
	yy <- (yy+min(data[[i]]$x))
	yy[yy < 0] <- 0
	if (slope >1){
		xx2 <- c(xx[1], xx, xx[length(xx)])
		yy2 <- c(0, yy, 0)
		if(plotAUC){
			polygon(xx2, yy2, density=15, col="red")
			}
		points(xx, yy, type="l", col="black", lwd=2)					
		if(ZOI == 20){
			points(xx20, yy20halo, col="blue4", cex=2, pch=19)
			}
		if(ZOI ==50){
			points(xx50, yy50halo, col="blue", cex=2, pch=19)
			}
		if(ZOI ==80){
			points(xx80, yy80halo, col="deepskyblue", cex=2, pch=19)
			}
		if(ZOI=="all"){
			points(xx80, yy80halo, col="blue4", cex=1.75, pch=19)
			points(xx50, yy50halo, col="blue", cex=1.75, pch=19)
			points(xx20, yy20halo, col="deepskyblue", cex=1.75, pch=19)
			}			
		}
	mtext(label, side=3, cex=0.6)
}

.plotAUC <- function(projectName, ML , ML2, stand,  clearHaloStand, standardLoc = 2.5, ymax=200, dotedge = 3.4, maxDist= 40, xplots = 4, height = 10, width=7,  AUC=50, ZOI=50, overwrite = TRUE, popUp = TRUE, plotAUC = TRUE, label=label, savePDF = TRUE, plotSub = plotSub){
	data <- eval(parse(text=projectName))
	if(is.na(plotSub[1])){
		plotSub <- 1:length(data)
		}
	fileFolder <- projectName
	dir.create(file.path(getwd(), "figures"), showWarnings= FALSE)
	dir.create(file.path(getwd(), "figures", fileFolder), showWarnings= FALSE)
	t <- file.path("figures", projectName , paste(projectName, "_AUC.pdf", sep=""))
	if (!overwrite){
		if (file.exists(t)){
			t <- file.path("figures", projectName , paste(projectName, "_AUC_2_AUC", AUC, "_ZOI", ZOI, ".pdf", sep=""))
			if (file.exists(t)){
				k <- 2
				while(file.exists(t)){
					k <- k+1
					t <- file.path("figures", projectName, paste(projectName, "_AUC_", k, "_AUC", AUC, "_ZOI", ZOI, ".pdf", sep=""))
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
			.singleAUC(data = data, ML = ML, ML2 = ML2, dotedge = dotedge, maxDist = maxDist, ymax = ymax, stand = stand, i = k,AUC=AUC, ZOI = ZOI, clearHaloStand = clearHaloStand, label=label[k], plotAUC = plotAUC)
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
