#' Add factor column
#' @description Add an extra factor ("type") column to the existing dataframe

#' @inheritParams maxLik
#' @param typePlace a number that indicates the position of the photograph name to be stored as the 'type' vector'. Defaults to 3. For more details see \code{\link{IJMacro}}
#' @param typeName a character string that indicates what to name the typeVector. Defaults to "type2".

#' @return updates the existing dataframe 

#' @export

#' @author Aleeza C. Gerstein


addType <- function(projectName, typePlace=3, typeName="type2"){
	df <- eval(parse(text=paste(projectName, ".df", sep="")))
	type2 <- unlist(lapply(df$name, function(x) strsplit(as.character(x), "_")[[1]][typePlace]) )
	place <- which.max(names(df) == "ZOI80")
	dfnew <- cbind(df[,1:place-1], type2, df[,place:length(df)])
	if(typeName != "type2") names(dfnew)[place] <- typeName

	filename <- file.path(getwd(), "parameter_files", paste(projectName, "_df.csv", sep="")	)
	dfName <- paste(projectName, ".df", sep="")
	cat(paste(dfName, " has been written to the global environment", sep=""))
	cat(paste("\nSaving file: ", filename,  sep=""))
	cat(paste("\n", projectName, "_df.csv can be opened in MS Excel.",  sep=""))
	assign(dfName, dfnew, envir=globalenv())
	}