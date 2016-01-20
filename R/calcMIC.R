#' Convert RAD to MIC based on built-in or provided parameters or datasets

#' @description Used to convert RAD into MIC. This conversion can be based on a) existing built-in data from a number of species/drug combinations, b) a user-supplied slope and intercept of the linear relationship between RAD and log2(MIC) for their species/drug combination of interest, c) a user supplied file containing MIC information from lines previously analyzed by diskImageR for RAD, or d) a user supplied file containing both RAD and MIC information. 

#' @param projectName the short name you have been using for the project.
#' @inheritParams twoParamPlot
#' @inheritParams oneParamPlot
#' @param RAD a numeric value the the critical level of the radius of inhibition (i.e., resistance) parameter to use for MIC. Currently only \code{RAD} = "80" (80\% reduction in growth), \code{RAD} = "50" (50\% reduction in growth), and \code{RAD} = "20" (20\% reduction in growth) are supported [Default = "20"]
#' @param addBreakpoints Indicates whether to add breakpoint lines to the standard curve plot (if the user has supplied data to generate a standard curve)

#' @return In all cases the function will return an updated .csv file that contains the MIC values that correspond to calculated RAD values in the directory "parameter_files" in the main project directory. If the user has supplied their own MIC data the function will also save the calculated model parameters into a separate file and will plot the linear relationship and line of best fit.

#' @export

#'@author Aleeza C. Gerstein & Inbal Hecht


calcMIC <- function(projectName, type="df", RAD="20", height = 4, width = 6, addBreakpoints = TRUE, savePDF = TRUE, popUp = TRUE){
	ZOIvalue <- RAD
	knownSppDrug <- data.frame(number = c(1:3), species=rep("C. albicans", 3), drug=c("fluconazole", "voriconazole", "posaconazole"), intercept=c(10, 11, 12), slope=c(-0.4, -0.5, -0.6))	
	# knownSppDrug <- file.path(.libPaths(), "diskImageR", "knownSppDrug.csv")	
	#Check whether file exists in environment or prompt user to load it
	if(type == "df"){
		if(paste(projectName, ".df", sep="") %in% ls(globalenv())) dataframe <- eval(parse(text=paste(projectName, ".df", sep="")))
		else stop(paste(projectName, ".df not found in working environment. Please load with function 'readExistingDF'", sep=""))
		}
	if(type == "ag"){
		if(paste(projectName, ".ag", sep="") %in% ls(globalenv())) dataframe <- eval(parse(text=paste(projectName, ".ag", sep="")))
		else stop(paste(projectName, ".ag not found in working environment. Please load with function 'readExistingAG'", sep=""))	
		}	

	useBuiltIn <- readline("Do you want to use built-in data for existing species/drug combinations? [y/n]  ")
	if(useBuiltIn=="y"){
		#do this eventually
		print(knownSppDrug[,1:3])
		sppDrug <- readline("Please choose the number that corresponds to the species/drug combination you wish to use: ")
		curvePars <- c(knownSppDrug[sppDrug, 4], knownSppDrug[sppDrug, 5])
		cat(paste("\nintercept: ", curvePars[1], "\tslope: ", curvePars[2], "\n", sep=""))		
	}
	else{
		haveSI<-readline("Do you know the slope and intercept of the linear relationship between RAD and log2(MIC)? [y/n]  ")
		if(haveSI == "y"){
	   		intercept <- readline("Enter the intercept: ")
	   		intercept <- as.numeric(intercept)
	   		slope <- readline("Enter the slope: ")
	   		slope <- as.numeric(slope)
	   		curvePars <- c(intercept, slope)
			cat(paste("\nintercept: ", curvePars[1], "\tslope: ", curvePars[2], "\n", sep=""))		
	   	}
		if(haveSI =="n"){
	  		exFile <- readline("Do you have a standard curve file?  [y/n] ")   
	  		if(exFile == "n") stop("You must use data from built-in combinations, provide the required parameters or a standard curve file to proceed.")
		  	else{
 		 		sameData <- readline("Does your file have MIC data that corresponds to the same strains as your current dataset? [y/n] ")
				if(sameData == "y"){
					MICFile <- tcltk::tk_choose.files(caption = "Select the MIC standard curve file (text file, comma delimited)") 
					MICdata<- read.csv(MICFile, header=TRUE,sep=",") 
					if (ncol(MICdata)<2) stop("Wrong data format: the file must contain two columns, the first containing the line name, and the second wtih corresponding MIC values \n")
					else{
						MIC_names <- MICdata[,1]
						MIC<-0
						N1 <- dataframe[,1]
						for (i in 1:length(N1)){  
							ind <- grep(N1[i],MIC_names);
							if (length(ind)>0) MIC[i]<-mean(MICdata[ind,2])
							else MIC[i]<-NA
						}
						RAD <- subset(dataframe, names %in% MIC_names)[paste("RAD", RAD, sep="")]
						if(length(RAD) != length(MIC) & type=="df") stop(paste("The length of the MIC file does not match the length of ", projectName, ".df. If you have replicates for RAD please run 'aggregateData' and the rerun 'calcMIC' with 'type=\"ag\"'. If you have replicates for MIC please average in the standard curve file before proceeding", sep=""))
						if(length(RAD) != length(MIC) & type=="df") stop(paste("The length of the MIC file does not match the length of ", projectName, ".ag.", sep=""))						
					}
					fit<-lm(log(MIC)~RAD, na.action=na.exclude)
					A <- summary(fit)$coefficients[1]
					B <- summary(fit)$coefficients[2]
					curvePars <-c(A, B)
					cat(paste("\nintercept: ", round(curvePars[1],2), "\tslope: ", round(curvePars[2],2), "\n", sep=""))
				}
				if(sameData == "n"){
					MICFile <- tcltk::tk_choose.files(caption = "Select the MIC standard curve file (text file, comma delimited)") 		
					MICdata <- read.csv(MICFile, header=TRUE)
					if (ncol(MICdata)<2) stop("Wrong data format: the file must contain at least two columns, one containing RAD values ('RAD'), and one with corresponding MIC values ('MIC') \n")
					else{
						RAD <- MICdata$RAD
						MIC <- MICdata$MIC
						fit <- lm(log2(MIC)~RAD, na.action=na.exclude)
						A <- summary(fit)$coefficients[1]
						B <- summary(fit)$coefficients[2]
						curvePars <-c(A, B)
						cat(paste("\nintercept: ", round(curvePars[1],2), "\tslope: ", round(curvePars[2],2), "\n", sep=""))		
						}
					}
				}			
		if(savePDF){
			t <- file.path(getwd(), "figures", projectName,  "RAD-MIC_standardCurve.pdf")					
			pdf(t, width=width, height=height)
			par(oma=c(1, 4, 1, 1))
			plot(RAD, log2(MIC), xlim=c(0, max(RAD)+3), yaxt="n", xaxt="n", xlab="", ylab="")
			axis(1, cex.axis=0.8)
			axis(2, at=log2(c(0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)), las=2, labels=c("0.12", "0.25", "0.5", "1", "2", "4", "8", "16", "32", "64", "128"), cex.axis=0.8)
			if(addBreakpoints){
				abline(v=(19-6)/2, col="grey")
				abline(v=(14.5-6)/2, col="grey")
				abline(h=log2((32+64)/2), col="grey")
				abline(h=log2((8+16)/2), col="grey")
				}
			abline(lm(log2(MIC)~RAD), col="red")
			mtext("Disk zone (mm)", side=1, outer=FALSE, line=3) 			
			title <- as.list(expression(paste(log[2], "(MIC)", sep=""), paste("actual values indicated, ", mu, "g/mL)")))
			mtext(do.call(expression, title), side=2, cex=0.8, line = c(3.5, 2.5))		
			dev.off()
			cat(paste("\nFigure saved: ", t, sep=""))
			if(popUp){
				tt <- paste("open ",t)
				system(tt)
			}
	}	
			paramName <- file.path(getwd(), "parameter_files", projectName, paste("RAD-MIC_parameters.csv", sep=""))
			params <- data.frame(intercept = curvePars[1], slope = curvePars[2])
			write.csv(params, file=paramName, row.names=FALSE)	
			cat(paste("\n\nThe calculated parameters have been saved here: \n", paramName, "\n and can be used in the future for the same species/drug combination\n", sep=""))
			}
		}
	MIC <- c(round(2^curvePars[1]*2^(curvePars[2]*dataframe[paste("RAD", ZOIvalue, sep="")]), 2)	)
	upDataframe <- data.frame(dataframe, MIC = MIC)
	dfName <- paste(projectName, ".df", sep="")	
	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df.csv", sep=""))
	write.csv(upDataframe, file=filename, row.names=FALSE)	
	cat(paste("\n", dfName, " has been updated and written to the global environment", sep=""))
	cat(paste("\n\nSaving file: ", filename, sep=""))
	assign(dfName, upDataframe, envir=globalenv())
	}