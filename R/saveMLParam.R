#' Save maximum likelihood output

#' @description Saves the output of maximum likelihood functions - asym, od50, scal, sigma, lnLik and type.

#' @inheritParams maxLik

#' @return A dataframe "projectName_ML.df" and a dataframe "projectName_ML2.df"are saved to the global environment and two .csv files "projectName_ML.csv" and "projectName_ML2.csv" are exported to the "parameter_files" directory when typical=TRUE, otherwise this function returns a dataframe "projectName_ML2.df" in global environment and a .csv file "projectName_ML2.csv" in "parameter_files" directory. 

#' @export

#' @author Aleeza C. Gerstein

saveMLparam <- function(projectName,typical=TRUE){
  fileFolder <- paste(Sys.Date(), projectName, sep="_")
  newdir <- file.path(getwd(), "parameter_files")
  newdir2 <- file.path(getwd(), "parameter_files", projectName)
  if (!file.exists(newdir)){		
    dir.create(newdir, showWarnings = FALSE)
    cat(paste("\nCreating new directory: ", newdir), sep="")
  }
  if (!file.exists(newdir2)){		
    dir.create(newdir2, showWarnings = FALSE)
    cat(paste("\nCreating new directory: ", newdir2), sep="")
  }
  
  if(typical){
    ML.df <- .MLparam(projectName)
    MLdf <- paste(projectName, "_ML.df", sep="")
    filename1 <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML.csv", sep=""))
    cat("\n")
    cat(paste("\n", MLdf, " has been written to the global environment", sep=""))
    assign(MLdf, ML.df, inherits=TRUE)
    cat(paste("\nSaving file: ", filename1))
    write.csv(ML.df, file=filename1, row.names=FALSE)
  }
  
  ML2.df <- .ML2param(projectName)
  ML2df <- paste(projectName, "_ML2.df", sep="")
  
  filename2 <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML2.csv", sep=""))	
  
  
  cat(paste("\n", ML2df, " has been written to the global environment", sep=""))
  assign(ML2df, ML2.df, inherits=TRUE)
  cat(paste("\nSaving file: ", filename2))
  write.csv(ML2.df, file=filename2, row.names=FALSE)	
}

.MLparam <- function(projectName){
  data <- eval(parse(text=projectName))
  ML <- eval(parse(text=paste(projectName, ".ML", sep="")))
  #unlist() simplifies the list into a single vector
  asym <- round(unlist(lapply(ML, function(x) x$par[1])), 2)
  od50 <- round(unlist(lapply(ML, function(x) x$par[2])), 2)
  scal <- round(unlist(lapply(ML, function(x) x$par[3])), 2)
  sigma <- round(unlist(lapply(ML, function(x) x$par[4])), 2)
  lnLik <- round(unlist(lapply(ML, function(x) x$lnLik)), 2)
  #line is the first column of data frame
  ML.df <- data.frame(line = names(data), asym, od50, scal, sigma, lnLik)
  return(ML.df)
}
#when typical=true, type column will be NA amounts.
.ML2param <- function(projectName){
  data <- eval(parse(text=projectName))
  ML2 <- eval(parse(text=paste(projectName, ".ML2", sep="")))
  #lapply: extracts the parametrs of ML2 list
  asymA <- round(unlist(lapply(ML2, function(x) x$par[1])), 2)
  od50A <- round(unlist(lapply(ML2, function(x) x$par[2])), 2)
  scalA <- round(unlist(lapply(ML2, function(x) x$par[3])), 2)
  sigma <- round(unlist(lapply(ML2, function(x) x$par[4])), 2)
  asymB <- round(unlist(lapply(ML2, function(x) x$par[5])), 2)
  od50B <- round(unlist(lapply(ML2, function(x) x$par[6])), 2)
  scalB <- round(unlist(lapply(ML2, function(x) x$par[7])), 2)
  lnLik <- round(unlist(lapply(ML2, function(x) x$lnLik)), 2)
  type <- unlist(lapply(ML2, function(x) x[[9]][["type"]]))
  ML2.df <- data.frame(line = names(data), asymA, od50A, scalA, sigma, asymB, od50B, scalB, lnLik,type)
  return(ML2.df)
}

