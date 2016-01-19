#' Convert RAD to MIC based on built-in or provided parameters or datasets

#' @description Used to convert RAD into MIC based on a standard curve.

#' @param projectName the short name you have been using for the project.
#' @inheritParams twoParamPlot
#' @param RAD a numeric value the the critical level of the radius of inhibition (i.e., resistance) parameter to use for MIC. Currently only \code{RAD} = "80" (80\% reduction in growth), \code{RAD} = "50" (50\% reduction in growth), and \code{RAD} = "20" (20\% reduction in growth) are supported [Default = "20".

#' @return An updated .csv file that contains the MIC values that correspond to calculated RAD values is saved to the directory "parameter_files" in the main project directory. The average line for each photograph is saved to the list \code{projectName} in the global environment.

#' @export

#'@author Aleeza C. Gerstein


calcMIC <- function(projectName, type="df", RAD="20"){
	knownSppDrug <- data.frame(number = c(1:3), species=rep("C. albicans", 3), drug=c("fluconazole", "voriconazole", "casitone"), intercept=c(10, 11, 12), slope=c(-0.4, -0.5, -0.6))	
	#probably need to check whether this exists in the global environment first, otherwise prompt user to open it.
	if(type == "df")
		if(paste(projectName, ".df", sep="") %in% ls()) data <- eval(parse(text=paste(projectName, ".df", sep="")))
		else stop(paste(data, ".df not found in working environment. Please load with function 'readExistingDF'", sep=""))
	if(type == "ag"){
		if(paste(projectName, ".ag", sep="") %in% ls()) data <- eval(parse(text=paste(projectName, ".ag", sep="")))
		else stop(paste(data, ".ag not found in working environment. Please load with function 'readExistingAG'", sep=""))	
		}	

	dataframe <- eval(parse(text=paste(projectName, ".df", sep="")))
	useBuiltIn <- readline("Do you want to use built-in data for existing species/drug combinations? (y/n)  ")
	if(useBuiltIn=="y"){
		knownSppDrug <- file.path(.libPaths(), "diskImageR", "knownSppDrug.csv")
		print(knownSppDrug[,1:3])
		sppDrug <- readline("Please choose the number that corresponds to the species/drug combination you wish to use: ")
		curvePars <- c(knownSppDrug[sppDrug, 4], knownSppDrug[sppDrug, 5])
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
	  		exFile <- readline("Do you have a standard curve file?  (y/n)\n ")   
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
						fit<-lm(log(MIC)~R1, na.action=na.exclude)
						A <- summary(fit)$coefficients[1]
						B <- summary(fit)$coefficients[2]
						curvePars <-c(exp(A), B)
					}
				}
				if(sameData == "n"){
					MICFile <- tcltk::tk_choose.files(caption = "Select the MIC standard curve file (text file, comma delimited)") 				
					if (ncol(MICdata)<2) stop("Wrong data format: the file must contain two columns, one containing RAD values ('RAD'), and one with corresponding MIC values ('MIC') \n")
					else{
						fit <- lm(log(MIC)~RAD, na.action=na.exclude)
						A <- summary(fit)$coefficients[1]
						B <- summary(fit)$coefficients[2]
						curvePars <-c(exp(A), B)
						}
					}
				}
			}
		}
			cat(paste(expression("The calculated equation from your standard curve is : ", log[2], "MIC = ", curvePars[1], "x", 2^(curePars[2]*RAD), sep="")))
			cat(" The parameters have been saved to XXXXX, which can be used in the future for the same species/drug combination")
			cat (" A figure of this standard curve has been saved to XXXX")  
	          #ADD in plot 		
	dataframe$MIC<-intercept*exp(slope*dataframe[paste("RAD", RAD, sep="")])
	dfName <- paste(projectName, ".df", sep="")	
	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df.csv", sep=""))
	write.csv(dataframe, file=filename, row.names=FALSE)	
	cat(paste("\n", dfName, " has been written to the global environment", sep=""))
	cat(paste("\n\nSaving file: ", filename, sep=""))
	assign(dfName, ag, envir=globalenv())
	}

calcMICown <- function(projectName, RAD="20"){
	dataframe <- eval(parse(text=paste(projectName, ".df", sep="")))
 
	  	
    	
	intercept <- knownSppDrug[sppDrug, 4]
	slope <- knownSppDrug[sppDrug, 5]
	dataframe$MIC<-intercept*exp(slope*dataframe[paste("RAD", RAD, sep="")])
	dfName <- paste(projectName, ".df", sep="")	
	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df.csv", sep=""))
	write.csv(dataframe, file=filename, row.names=FALSE)	
	cat(paste("\n", dfName, " has been written to the global environment", sep=""))
	cat(paste("\n\nSaving file: ", filename, sep=""))
	assign(dfName, ag, envir=globalenv())
	}


        
# +          if(exFile=="y"){ 
# +            setwd(getwd())
# +             ;
# +            tempPars<-read.table(curveFile, header=FALSE,sep=",",dec=".");
# +            CurvePars=unlist(tempPars)
# +
# +            
# +            # calculate MIC and add to data file
# +            }
# +          else{
# +            calib<-readline("Do you want to provide data for MIC calibration? (y/n)\n")
# +            if(calib=="y"){
# +               print("Please select the file in the opened browser",quote = FALSE)                
# +               MICFile <- tcltk::tk_choose.files(caption = "Select the MIC calibration file: Text file, tab delimited") 
# +               MICdata<- read.csv(MICFile, header=FALSE,sep="\t",dec="."); 
# +             +               MIC_length=length(MICdata[[1]])
# +               MIC_names=MICdata$V1
# +              # R1<-0
# +               MIC<-0
# +               for (i in 1:length(N1))
# +                 {
# +                 #R1[i]<-RAD2.df[MICdata[i,1],"RAD20"];
# +                 ind=grep(N1[i],MIC_names);
# +                 if (length(ind)>0)
# +                    MIC[i]<-mean(MICdata[ind,2])
# +                 else
# +                   MIC[i]<-NA
# +                  }
# +               fit<-lm(log(MIC)~R1,na.action=na.exclude)
# +               A=summary(fit)$coefficients[1]
# +               CurvePars<-exp(A)
# +               B=summary(fit)$coefficients[2]
# +               CurvePars[2]<-B
# +               CurveFile <- readline( "How would you like to name the curve file?\n ") 
# +               CurveFile=paste(getwd(),"/",CurveFile,".csv",sep="")
# +               # handle the case of a pre-existing file
# +               OV_flag=FALSE;
# +               while (!OV_flag && TRUE %in% file.exists(CurveFile))
# +               {
# +                 ov<-readline("File already exists, overwrite? (y/n) \n")
# +                 if (ov=="y") OV_flag=TRUE
# +                 else {
# +                   CurveFile <- readline( "How would you like to name the curve file?  \n ")
# +                   CurveFile=paste(getwd(),"/",CurveFile,".csv",sep="")
# +                 }
# +                       
# +               }
# +               write.table(CurvePars,CurveFile,col.names = FALSE,row.names = FALSE)
# +              }
# +            else {CurvePars[1]=0
# +                  CurvePars[2]=0}
# +          }
# +        }
# +    }
# +   CurvePars    
# +  }
# +
# + 
# +
# +
# diff --git a/R/createDataframe.R b/R/createDataframe.R
# index b99c822..b4a7809 100644
# --- a/R/createDataframe.R
# +++ b/R/createDataframe.R
# @@ -97,8 +97,18 @@ createDataframe <- function(projectName, clearHalo, diskDiam = 6, maxDist = 30,
 	# aveFoG80 <- FoG80/x80
 	# aveFoG50 <- FoG50/x50
 	# aveFoG20 <- FoG20/x20	
# -
# +  
 	# param <- data.frame(RAD80 =round(x80, digits=0), RAD50 = round(x50, digits=0), RAD20 = round(x20, digits=0), FoG80 = round(FoG80/maxFoG80, digits=2), FoG50 = round(FoG50/maxFoG50, digits=2), FoG20 = round(FoG20/maxFoG20, digits=2), slope=round(slope, digits=1))
 	
# +	#### Calculate MIC: function calc_curve imports or calculates curve parameters based on the user's choice and input
# +	exp_names <- unlist(names(data))
# +	curve=calc_curve(exp_names,x20)
# +	RAD20=round(x20,digits=0)
# +
# +	MIC=curve[1]*exp(curve[2]*RAD20)
# +	# MIC is added to the output Excel file only if the user wants to calculate it
# +	if (curve[1]>0)
# +	  param<-data.frame(param,MIC=round(MIC,digits=2))
# +	###
 	
 	# if (is.logical(nameVector)){
 		# if (nameVector){
