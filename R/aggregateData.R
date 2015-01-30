#' Averages photographs of the same type

#' @description Uses a user-supplied variance measure (currently supported: standard error, coefficient of variation, built-in R functions (e.g., sd) to calculate variance among photographs of the same type

#' @inheritParams maxLik
#' @param replicate a character vector indicating which the column names that contain which factors to use. Defaults to c("line", "type"). Note that if the typeVector name was changed in \code{createDataframe} this should be reflected here.
#' @param varFunc what type of variation measurment to perform. Currently supports \code{varFunc} = "se" to calculate the standard error, \code{varFun} = "cv" to calculate the coefficient of variation or any built-in R function (e.g., sd). 
#' @param overwrite a logical value indicating whether to overwrite existing aggregate dataframe for the same project name. This allows you to save different dataframes averaging across different factors or using different variance measures

#' @return A dataframe "projectName.ag" is saved to the global environment and a .csv file "projectName_ag.csv" is exported to the "parameter_files" directory. 

#' @export

#' @author Aleeza C. Gerstein


aggregateData <- function(projectName, varFunc = "se", replicate = c("line", "type"), overwrite = TRUE){
	dataframe <- eval(parse(text=paste(projectName, ".df", sep="")))
	
	if (varFunc == "se") var <- se
	if (varFunc == "cv") var <- cv
	if (!varFunc %in% c("se", "cv"))  var <- eval(parse(text=varFunc))
	
	temp <- aggregate(dataframe[c("ZOI80", "ZOI50", "ZOI20", "fAUC80", "fAUC50", "fAUC20", "slope")], dataframe[replicate], mean, na.rm=TRUE)
	t <- apply(dataframe[,c("ZOI80", "ZOI50", "ZOI20", "fAUC80", "fAUC50", "fAUC20", "slope")], 2, function(x) aggregate(x, dataframe[replicate], var))
	 var <- data.frame(t[[1]]$x, t[[2]]$x, t[[3]]$x, t[[4]]$x, t[[5]]$x, t[[6]]$x, t[[7]]$x)
	names(var) <- paste(varFunc, names(t), sep=".")
	ag <- cbind(temp, var)
	
	ag[,names(ag) %in% c("ZOI80", "ZOI50", "ZOI20", "slope")] <- round(ag[,names(ag) %in% c("ZOI80", "ZOI50", "ZOI20", "slope")])
	ag[,	names(ag) %in% c("fAUC80", "fAUC50", "fAUC20", paste(varFunc, "ZOI80", sep="."), paste(varFunc, "ZOI50", sep="."), paste(varFunc, "ZOI20", sep="."), paste(varFunc, "slope", sep="."))] <- round(ag[,names(ag) %in% c("fAUC80", "fAUC50", "fAUC20", paste(varFunc, "ZOI80", sep="."), paste(varFunc, "ZOI50", sep="."), paste(varFunc, "ZOI20", sep="."), paste(varFunc, "slope", sep="."))], digits=2)	
	ag[, names(ag) %in% c(paste(varFunc, "fAUC80", sep="."), paste(varFunc, "fAUC50", sep="."), paste(varFunc, "fAUC20", sep="."))] <- round(ag[, names(ag) %in% c(paste(varFunc, "fAUC80", sep="."), paste(varFunc, "fAUC50", sep="."), paste(varFunc, "fAUC20", sep="."))], digits=4)	
	
	newdir2 <- file.path(getwd(), "parameter_files", sep="")		
	newdir3 <- file.path(getwd(), "parameter_files", projectName)	
	dir.create(newdir2, showWarnings = FALSE)
	dir.create(newdir3, showWarnings = FALSE)
	agName <- paste(projectName, ".ag", sep="")
	filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ag.csv", sep=""))

	if (!overwrite){
		if (file.exists(filename)){
			filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ag_2.csv", sep=""))
			agName <- paste(projectName, "_2.ag", sep="")
			if (file.exists(filename)){
				k <- 2
				while(file.exists(filename)){
					k <- k+1
					filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ag_", k, ".csv", sep=""))
					agName <- paste(projectName, "_", k, ".ag", sep="")	
					}
				}
			}
		}

	
	write.csv(ag, file=filename, row.names=FALSE)	

	cat(paste("\n", agName, " has been written to the global environment", sep=""))
	cat(paste("\n\nSaving file: ", filename, sep=""))
	
	 assign(agName, ag, envir=globalenv())
	}
	
cv <- function(x, na.rm=TRUE) (100*sd(x)/mean(x))
se <- function(x, na.rm=TRUE) sd(x)/sqrt(length(x))