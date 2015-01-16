#' Averages photographs of the same type

#' @description Addin description here

#' @inheritParams maxLik
#' @param replicate a character vector indicating which the column names that contain which factors to use. Defaults to c("lines", "type"). Note that if the typeVector name was changed in \code{createDataframe} this should be reflected here.
#' @param varFunc what type of variation measurment to perform. Currently supports \code{varFunc} = "se" to calculate the standard error, \code{varFun} = "cv" to calculate the coefficient of variation or any built-in R function (e.g., sd). 

#' @return A dataframe "projectName.ag" is saved to the global environment and a .csv file "projectName_ag.csv" is exported to the "parameter_files" directory. 

#' @export

#' @author Aleeza C. Gerstein


aggregateData <- function(projectName, replicate = c("lines", "type"), varFunc = "se"){
	dataframe <- eval(parse(text=paste(projectName, ".df", sep="")))
	
	if (varFunc == "se") var <- .se

	if (varFunc == "cv") var <- .cv
	else var <- eval(parse(text=varFunc))
	temp <- aggregate(c(dataframe[4:10]), dataframe[replicate], mean, na.rm=TRUE)
	t <- apply(dataframe[,4:10], 2, function(x) aggregate(x, dataframe[replicate], var))
	var <- data.frame(t[[1]]$x, t[[2]]$x, t[[3]]$x, t[[4]]$x, t[[5]]$x, t[[6]]$x, t[[7]]$x)
	names(var) <- paste(varFunc, names(t), sep=".")
	ag <- cbind(temp, var)

	ag[,c(3:5, 9)] <- round(ag[,c(3:5, 9)])
	ag[,c(6:8, 10:12, 16)] <- round(ag[,c(6:8, 10:12, 16)], digits=2)	
	ag[,c(13:15)] <- round(ag[,c(13:15)], digits=4)	
	
	filename <- paste(getwd(), "/parameter_files/", Sys.Date(), "_", projectName, "/",projectName, "_ag.csv", sep="")
	newdir2 <- paste(getwd(), "/parameter_files/", sep="")		
	newdir3 <- paste(getwd(), "/parameter_files/", Sys.Date(), "_", projectName, "/", sep="")	

	dir.create(newdir2, showWarnings = FALSE)
	dir.create(newdir3, showWarnings = FALSE)

	write.csv(ag, file=filename, row.names=FALSE)	
	agName <- paste(projectName, ".ag", sep="")
	cat(paste("\n", agName, " has been written to the global environment", sep=""))
	cat(paste("\n\nSaving file: ", filename, sep=""))
	cat(paste("\n",  projectName, "_ag.csv can be opened in MS Excel (save as .xls file if desired)",  sep=""))
	assign(agName, ag, envir=globalenv())
	}

.cv <- function(x, na.rm=TRUE) (100*sd(x)/mean(x))
.se <- function(x, na.rm=TRUE) sd(x)/sqrt(length(x))