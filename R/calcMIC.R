#' Convert RAD to MIC based on built-in or provided parameters or datasets

#' @description Used to convert RAD into MIC. This conversion can be based on a) existing built-in data from a number of species/drug combinations, b) a user-supplied slope and intercept of the linear relationship between RAD and log2(MIC) for their species/drug combination of interest, c) a user supplied file containing MIC information from lines previously analyzed by diskImageR for RAD, or d) a user supplied file containing both RAD and MIC information. 

#' @param projectName the short name you have been using for the project.
#' @inheritParams twoParamPlot
#' @param RAD a numeric value the the critical level of the radius of inhibition (i.e., resistance) parameter to use for MIC. Currently only \code{RAD} = "80" (80\% reduction in growth), \code{RAD} = "50" (50\% reduction in growth), and \code{RAD} = "20" (20\% reduction in growth) are supported [Default = "20".
#' @param addBreakpoints Indicates whether to add breakpoint lines to the standard curve plot (if the user has supplied data to generate a standard curve).

#' @return In all cases the function will return an updated .csv file that contains the MIC values that correspond to calculated RAD values in the directory "parameter_files" in the main project directory. If the user has supplied their own MIC data the function will also save the calculated model parameters into a separate file and will plot the linear relationship and line of best fit.

#' @export

#'@author Aleeza C. Gerstein


calcMIC <- function(projectName, type="df", RAD="20", addBreakpoints = TRUE, savePDF = TRUE, popUp = TRUE){
	knownSppDrug <- data.frame(number = c(1:3), species=rep("C. albicans", 3), drug=c("fluconazole", "voriconazole", "casitone"), intercept=c(10, 11, 12), slope=c(-0.4, -0.5, -0.6))	
	#Check whether file exists in environment or prompt user to load it
	if(type == "df"){
		if(paste(projectName, ".df", sep="") %in% ls(globalenv())) dataframe <- eval(parse(text=paste(projectName, ".df", sep="")))
		else stop(paste(projectName, ".df not found in working environment. Please load with function 'readExistingDF'", sep=""))
		}
	if(type == "ag"){
		if(paste(projectName, ".ag", sep="") %in% ls(globalenv())) dataframe <- eval(parse(text=paste(projectName, ".ag", sep="")))
		else stop(paste(projectName, ".ag not found in working environment. Please load with function 'readExistingAG'", sep=""))	
		}	

	useBuiltIn <- readline("Do you want to use built-in data for existing species/drug combinations? (y/n)  ")
	if(useBuiltIn=="y"){
		#do this eventually
		# knownSppDrug <- file.path(.libPaths(), "diskImageR", "knownSppDrug.csv")
		print(knownSppDrug[,1:3])
		sppDrug <- readline("Please choose the number that corresponds to the species/drug combination you wish to use: ")
		curvePars <- c(knownSppDrug[sppDrug, 4], knownSppDrug[sppDrug, 5])
		# cat(paste(" interecept:," curvePars[1], "slope:", curvePars[2], sep=""))
	}
	else{
		haveSI<-readline("Do you know the slope and intercept of the linear relationship between RAD and log2(MIC)? (y/n)  ")
		if(haveSI == "y"){
	   		intercept <- readline("Enter the intercept: ")
	   		intercept <- 2^intercept
	   		slope <- readline("Enter the slope: ")
	   		curvePars <- c(intercept, slope)
	   	}
		if(haveSI =="n"){
	  		exFile <- readline("Do you have a standard curve file?  (y/n) ")   
	  		if(exFile == "n") stop("You must use data from built-in combinations, provide the required parameters or a standard curve file to proceed.")
		  	else{
 		 		sameData <- readline("Does your file have MIC data that corresponds to the same strains as your current dataset? (y/n)")
				if(sameData == "y"){
					MICFile <- tcltk::tk_choose.files(caption = "Select the MIC standard curve file (text file, comma delimited)") 
					MICdata<- read.csv(MICFile, header=FALSE,sep=",") 
					if (ncol(MICdata)<2) stop("Wrong data format: the file must contain two columns, one containing the line name, and one wtih corresponding MIC values \n")
					else{
						MICdata <- read.csv(MICFile, header=FALSE,sep="\t",dec="."); 
						MIC_length <- length(MICdata[[1]])
						MIC_names <- MICdata$V1
						MIC<-0
						N1 <- dataframe[,1]
						for (i in 1:length(N1)){  
							ind <- grep(N1[i],MIC_names);
							if (length(ind)>0) MIC[i]<-mean(MICdata[ind,2])
							else MIC[i]<-NA
						}
						RAD <- subset(dataframe, names %in% MIC_names)[paste("RAD", RAD, sep="")]
						if(length(RAD) != length(MIC) & type=="df"){
							stop(paste("The length of the MIC file does not match the length of ", projectName, ".df. If you have replicates for RAD please run 'aggregateData' and the rerun 'calcMIC' with 'type=\"ag\"'. If you have replicates for MIC please average in the standard curve file before proceeding", sep=""))
						if(length(RAD) != length(MIC) & type=="df"){
							stop(paste("The length of the MIC file does not match the length of ", projectName, ".ag.", sep="")
						} 
						fit<-lm(log(MIC)~RAD, na.action=na.exclude)
						A <- summary(fit)$coefficients[1]
						B <- summary(fit)$coefficients[2]
						curvePars <-c(exp(A), B)
					}
				}
				if(sameData == "n"){
					MICFile <- tcltk::tk_choose.files(caption = "Select the MIC standard curve file (text file, comma delimited)") 				
					if (ncol(MICdata)<2) stop("Wrong data format: the file must contain two columns, one containing RAD values ('RAD'), and one with corresponding MIC values ('MIC') \n")
					else{
						RAD <- MICFile$RAD
						MIC <- MICFile$MIC
						fit <- lm(log(MIC)~RAD, na.action=na.exclude)
						A <- summary(fit)$coefficients[1]
						B <- summary(fit)$coefficients[2]
						curvePars <-c(exp(A), B)
						}
					}
				}			
			cat(paste(expression("\n The calculated equation from your standard curve is : ", log[2], "MIC = ", curvePars[1], "x", 2^(curePars[2]*RAD), sep="")))
	if(savePDF){
		pdf(t, width=width, height=height)
		plot(RAD, log2(MIC), xlim=c(0, 35), yaxt="n", xaxt="n", xlab="", ylab="")
		axis(1)
		axis(2, at=log2(unique(MIC)), las=2, labels=unique(MIC))
		if(addBreakpoints){
			abline(v=(19-6)/2)
			abline(v=(14.5-6)/2)
			abline(h=log2((32+64)/2))
			abline(h=log2((8+16)/2))
			}
		abline(lm(log2(MIC)~RAD), col="red")
		mtext("Disk zone (mm)", side=1, outer=FALSE, line=3) 			
		mtext("log2(MIC)\n (actual values indicated, ug/mL)", side=2, outer=FALSE, line=3.5) 			
		dev.off()
		t <- file.path(paste("figures", projectName,  "RAD-MIC_standardCurve.pdf", sep=""))		
		cat (paste("\nA figure of this standard curve has been saved to ", t, sep=""))  
		if(popUp){
			tt <- paste("open ",t)
			system(tt)
		}
	}
	
	        #ADD in plot
			paramName <- file.path(getwd(), "parameter_files", projectName, paste("RAD-MIC_parameters.csv", sep=""))
			params <- data.frame(intercept = curvePars[1], slope = curvePars[2])
			write.csv(params, file=paramName, row.names=FALSE)	
			cat(paste("\n\The calculated parameters have been saved here: ", paramName, "\n and can be used in the future for the same species/drug combination", sep=""))
			
			}
		}
	dataframe["MIC"] <- round(curvePars[1]*exp(curvePars[2]*dataframe[paste("RAD", RAD, sep="")]), 2)
	print(dataframe)
	dfName <- paste(projectName, ".df", sep="")	
	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df.csv", sep=""))
	write.csv(dataframe, file=filename, row.names=FALSE)	
	cat(paste("\n", dfName, " has been written to the global environment", sep=""))
	cat(paste("\n\nSaving file: ", filename, sep=""))
	assign(dfName, dataframe, envir=globalenv())
	}