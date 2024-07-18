 #' Dataframe creation

#' @description \code{createDataframe} saves the calculated resistance, perseverance and sensitivity estimates

#' @inheritParams maxLik
#' @param nameVector either a logical value or a character vector. Supported values are \code{nameVector} = "TRUE" to assign the photograph name to the 'name' column, \code{nameVector} = "FALSE" to assign the photograph number to the 'name' column, or \code{nameVector} = a vector the same length as the number of photographs indicating the desired names.
#' @param typeVector a logical value. \code{typeVector} = "TRUE" will add a 'type' vector to the dataframe using values found in the \code{typePlace} position of the photograph names (see \code{\link{IJMacro}} for more details) while \code{typeVector} = "FALSE" will not add a type column.
#' @param typePlace a number that indicates the position of the photograph name to be stored as the 'type' vector'. Defaults to 2. For more details see \code{\link{IJMacro}}
#' @param typeName a character string that indicates what to name the typeVector. Defaults to "type".
#' @param removeClear a logical value that indicates whether to remove the clear halo picture from the dataset (i.e., is this picture an experimental picture, or one solely included to use as a clear halo). Defaults to FALSE.
#' @param standType either "one" or "indiv" to determine whether to use one standard for all photos or individually standardize each photo.
#' @param standardLoc is a numberic value that indicates the location (on the disk) to use to standardize white intensity across photographs. The position of standardLoc is a position that should theoretically have the same intensity in all photographs, i.e., the white of the disk. The default value (2.5mm) was chosen after testing of 6mm disks that contain some writing. If smaller disks are used standardLoc should be scaled appropriately. You can see where standardLoc falls in each photograph in \code{plotRaw} (the red dashed line when `plotStandardLoc = TRUE`). To suppress this standardization use `standardLoc = FALSE`
#' @param typical if TRUE, a logistic curve will be calculated for each photo. If FALSE, the function will determine whether logistic, confounding, or paradoxical is the best fit and then calculate the curve.

#' @details If typical = TRUE, creates a dataframe with 9 columns: \bold{name:} determined by \code{nameVector}, either photograph names, photograph numbers, or a user-supplied list of names; \bold{line:} the first components of the \code{namesVector}, everything that comes before the first "_" in the photograph name; \bold{type:} the location within the \code{name} of the photograph type is supplied by \code{typePlace}. Use \code{\link{addType}} if more than one type column are desired; \bold{RAD20, RAD50, RAD80:} resistance parameters, corresponding to the distance in mm of 80\%, 50\% and 20\% reduction in growth; \bold{FoG80, FoG50, FoG20:} perseverance parameters, corresponding to the fraction of growth achieved above the 80\%, 50\% and 20\% reduction in growth points
#' 		
#' If typical = FALSE, then up to three dataframes will be created, depending on the categorizations of the pictures in maxLik. The pictures will be split into standard, confounding, and paradoxical growth and each category that has photos will create a dataframe. All the dataframes will have the same first three columns as typical = TRUE, as well as a column indicating the photo index so that photos can be identified. The parameters of interest vary between the drug responses.
#' The standard dataframes will have the same parameters of interest as typical = TRUE. 
#' The confounding dataframes will have DRAD80, DRAD50, and DRAD20, corresponding to the distance in mm of 80\%, 50\%, and 20\% disinhibition.
#' The paradoxical dataframes will have two points of maximum inhibition. One is the calculated maxInhib (based on the minimum point from the fit curve) and the other is the observed maxInhib (based on the minimum observed data point) 

#' @return If typical = TRUE, a dataframe "projectName.df" is saved to the global environment and a .csv file "projectName_df.csv" is exported to the "parameter_files" directory.
#' If typical = FALSE, up to three dataframes are saved to the global environment and up to three .csv files are exported to the "parameter_files" directory.

#' @examples
#' \dontrun{
#' createDataframe("myProject", clearHalo=1)
#' createDataframe("myProject", clearHalo=1, removeClear = TRUE, typeName = "drugAmt")
#' }

#' @export

createDataframe <- function(projectName, clearHalo, diskDiam = 6, maxDist = 25, standardLoc = 2.5, removeClear = FALSE, nameVector=TRUE, typeVector=TRUE, typePlace=2, typeName = "type", standType = "one", typical = TRUE){
  
  if(standType=="one"){
    if(!(hasArg(clearHalo))){
      cont <- readline(paste("Please specify photograph number with a clear halo ", sep=""))
      clearHalo <- as.numeric(cont)
    }
  }
  
  # standType = 'indiv' is not yet supported when typical = FALSE
  # this method finds the minimum pixel intensity in the area near the drug disk and uses it to standardize
  # the area near the drug disk often has heavy growth for confounding and paradoxical curves so this method can't be used
  if (!typical & standType == 'indiv') {
    stop('Individual standardization with paradoxical curves is not currently supported')
  }
  
  data <- eval(parse(text=projectName))
  
  df <- data.frame()
  dotedge <- diskDiam/2 + 0.7
  newdir <- file.path(getwd(), "parameter_files")
  newdir2 <- file.path(getwd(), "parameter_files", projectName)
  newdir3 <- file.path(getwd(), "figures", projectName)
  
  if (!file.exists(newdir)){
    dir.create(newdir, showWarnings = FALSE)
    cat(paste("\n\tCreating new directory: ", newdir), sep="")
  }
  if (!file.exists(newdir2)){
    dir.create(newdir2, showWarnings = FALSE)
    cat(paste("\nCreating new directory: ", newdir2), sep="")
  }
  if (!file.exists(newdir3)){
    dir.create(newdir3, showWarnings = FALSE)
    cat(paste("\nCreating new directory: ", newdir3), sep="")
  }
  
  # ML estimates will only be present if maxLik used typical = TRUE
  if (typical) {
    ML <- paste(projectName, ".ML", sep="")
    ML <- eval(parse(text=ML))
  }
  
  # if typical = FALSE, ML2 will contain all the parameter estimates
  ML2 <- paste(projectName, ".ML2", sep="")
  ML2 <- eval(parse(text=ML2))
  
  if (typical) {
    
    dotMax <- max(sapply(data, function(x) {x[which(x[,1] > standardLoc)[1], 2]}))
    stand <-c( sapply(data, function(x) {dotMax-x[which(x[,1] > standardLoc)[1], 2]}))
    
    if (standType == 'one') {
      
      clearHaloData <- data[[clearHalo]]
      # add the 0.5 to dotedge to be consistant with maxLik
      startX <- which(clearHaloData[,1] > dotedge+0.5)[1]
      stopX <- which(clearHaloData[,1] > maxDist - 0.5)[1]
      clearHaloData <- clearHaloData[startX:stopX, 1:2]
      clearHaloData$x <- clearHaloData$x + stand[clearHalo]
      clearHaloData$distance <- clearHaloData$distance - (dotedge+0.5)
      clearHaloStand <- clearHaloData[1,2]
      
      # .findFoG creates a list with all the calculated RAD/FoG values for each photo
      FoG.df <-  sapply(c(1:length(data)), .findFoG, data=data, ML=ML, ML2 = ML2, stand = stand, dotedge = dotedge,  maxDist = maxDist, clearHaloStand = clearHaloStand, standardLoc = standardLoc, standType = 'one')
      
    } else {
      
      clearHaloStand <- NULL
      FoG.df <-  sapply(c(1:length(data)), .findFoG, data=data, ML=ML, ML2 = ML2, stand = stand, dotedge = dotedge,  maxDist = maxDist, clearHaloStand = clearHaloStand, standardLoc = standardLoc, standType = 'indiv')
      
    }
    
    # need to unlist the parameters so that they can be put together nicely in a dataframe
    x80 <- unlist(FoG.df[1,])
    x50 <- unlist(FoG.df[2,])
    x20 <- unlist(FoG.df[3,])
    FoG80 <- unlist(FoG.df[4,])
    FoG50 <- unlist(FoG.df[5,])
    FoG20 <- unlist(FoG.df[6,])
    
    param <- data.frame(RAD80 = x80, RAD50 = x50, RAD20 = x20, FoG80 = FoG80, FoG50 = FoG50, FoG20 = FoG20)
    
    # addNames will add name/type/line columns to the dataframe 
    # this step was not done in a helper function in the legacy version
    param <- addNames(data, param, nameVector, typeVector, typePlace, typeName)
    
    # adjusting for edge cases
    # since FoG is the fraction of growth, it cannot be greater than 1
    param$FoG80[param$FoG80 >1] <- 1
    param$FoG50[param$FoG50 >1] <- 1
    param$FoG20[param$FoG20 >1] <- 1
    # if the RAD distance is 0, the FoG will be negligible
    param$FoG80[param$RAD80 == 0] <- NA
    param$FoG50[param$RAD50 == 0] <- NA
    param$FoG20[param$RAD20 == 0] <- NA
    param$FoG80[param$RAD80 == 1] <- NA
    param$FoG50[param$RAD50 == 1] <- NA
    param$FoG20[param$RAD20 == 1] <- NA
    
    if(standType == "one"){
      if (removeClear)	param <- param[-clearHalo,]
    }
    
    filename <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df.csv", sep=""))
    write.csv(param, file=filename, row.names=FALSE)
    
    dfName <- paste(projectName, ".df", sep="")
    cat(paste("\n", dfName, " has been written to the global environment", sep=""))
    cat(paste("\nSaving file: ", filename,  sep=""))
    cat(paste("\n", projectName, "_df.csv can be opened in MS Excel.",  sep=""))
    assign(dfName, param, inherits=TRUE)
    
  } 
  
  if (!typical) {
    print("Creating atypical dataframes")
    dotMax <- max(sapply(data, function(x) {x[which(x[,1] > standardLoc)[1], 2]}))
    stand <-c( sapply(data, function(x) {dotMax-x[which(x[,1] > standardLoc)[1], 2]}))
    
    # only standType = 'one' is supported here
    clearHaloData <- data[[clearHalo]]
    startX <- which(clearHaloData[,1] > dotedge+0.5)[1]
    stopX <- which(clearHaloData[,1] > maxDist - 0.5)[1]
    clearHaloData <- clearHaloData[startX:stopX, 1:2]
    clearHaloData$x <- clearHaloData$x + stand[clearHalo]
    clearHaloData$distance <- clearHaloData$distance - (dotedge+0.5)
    clearHaloStand <- clearHaloData[1,2]
    
    logIndices <- c()
    negLogIndices <- c()
    paraIndices <- c()
    
    # rather than creating subsets of data, we create vectors with the indices of photos of each type
    for (i in 1:length(ML2)) {
      
      if (ML2[[i]][length(ML2[[i]])] == 'Log') {
        logIndices <- c(logIndices, i)
        
      } else if (ML2[[i]][length(ML2[[i]])] == 'NegLog') {
        negLogIndices <- c(negLogIndices, i)
        
      } else {
        paraIndices <- c(paraIndices, i)
      }
    }
    
    # we go through a similar process to when typical = TRUE
    # we just do it more than once and with different helper functions/subsets of the data
    
    if (length(negLogIndices) != 0) {
      
      # only supply a subset of the full list of data to the .findNegLogParams function
      paramNegLog <- sapply(c(1:length(negLogIndices)), .findNegLogParams, data[negLogIndices], ML2[negLogIndices], stand[negLogIndices], clearHaloStand, dotedge = dotedge, maxDist = maxDist)
      
      x80 <- unlist(paramNegLog[1, ])
      x50 <- unlist(paramNegLog[2, ])
      x20 <- unlist(paramNegLog[3, ])
      
      paramNegLog <- data.frame(photo.index = negLogIndices, DRAD80 = x80, DRAD50 = x50, DRAD20 = x20)
      paramNegLog <- addNames(data[negLogIndices], paramNegLog, nameVector, typeVector, typePlace, typeName)
      
      dfName <- paste(projectName, "_confound.df", sep="")
      assign(dfName, paramNegLog, inherits=TRUE)
      cat(paste("\n", dfName, " has been written to the global environment", sep=""))
      
      filenameNegLog <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df_negLog.csv", sep=""))
      write.csv(paramNegLog, file = filenameNegLog, row.names = FALSE)
      
      cat(paste("\nSaving file: ", filenameNegLog,  sep=""))
      cat(paste("\n", projectName, "_confound_df.csv can be opened in MS Excel.\n",  sep=""))
   }
    
    if (length(logIndices) != 0) {
      
      # typical = TRUE and typical = FALSE use the same helper function to get parameters
      # the only difference is that typical = TRUE uses ML parameters, which typical = FALSE does not have
      # so ML = NULL is supplied to .findFoG in this case, which .findFoG can account for
      paramLog <- sapply(c(1:length(logIndices)), .findFoG, data[logIndices], NULL, ML2[logIndices], stand[logIndices], clearHaloStand, dotedge, maxDist, standardLoc, standType)
      
      x80 <- unlist(paramLog[1, ])
      x50 <- unlist(paramLog[2, ])
      x20 <- unlist(paramLog[3, ])
      FoG80 <- unlist(paramLog[4, ])
      FoG50 <- unlist(paramLog[5, ])
      FoG20 <- unlist(paramLog[6, ])
      
      paramLog <- data.frame(photo.index = logIndices, RAD80 = x80, RAD50 = x50, RAD20 = x20, FoG80, FoG50, FoG20)
      paramLog <- addNames(data[logIndices], paramLog, nameVector, typeVector, typePlace, typeName)
      
      paramLog <- paramLog[order(paramLog$line),]
      # make the same adjustments here to make the data make sense
      paramLog$FoG80[paramLog$FoG80 >1] <- 1
      paramLog$FoG50[paramLog$FoG50 >1] <- 1
      paramLog$FoG20[paramLog$FoG20 >1] <- 1
      paramLog$FoG80[paramLog$RAD80 == 0] <- NA
      paramLog$FoG50[paramLog$RAD50 == 0] <- NA
      paramLog$FoG20[paramLog$RAD20 == 0] <- NA
      paramLog$FoG80[paramLog$RAD80 == 1] <- NA
      paramLog$FoG50[paramLog$RAD50 == 1] <- NA
      paramLog$FoG20[paramLog$RAD20 == 1] <- NA
      
      dfName <- paste(projectName, ".df", sep="")
      assign(dfName, paramLog, inherits=TRUE)
      cat(paste("\n", dfName, " has been written to the global environment", sep=""))
      
      filenameLog <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df_log.csv", sep=""))
      write.csv(paramLog, file=filenameLog, row.names=FALSE)
      
      cat(paste("\nSaving file: ", filenameLog,  sep=""))
      cat(paste("\n", projectName, "_df.csv can be opened in MS Excel.\n",  sep=""))
    }
    
    if (length(paraIndices) != 0) {
      
      paramPara <- sapply(c(1:length(paraIndices)), .findParadoxicalParams, data[paraIndices], ML2[paraIndices], stand[paraIndices], clearHaloStand, dotedge = dotedge, maxDist = maxDist)
      
      xmin <- unlist(paramPara[1, ])
      xmin2 <- unlist(paramPara[2, ])
      
      paramPara <- data.frame(photo.index = paraIndices, curve.maxInhib = xmin, observed.maxInhib = xmin2)
      paramPara <- addNames(data[paraIndices], paramPara, nameVector, typeVector, typePlace, typeName)
      
      dfName <- paste(projectName, "_para.df", sep="")
      assign(dfName, paramPara, inherits=TRUE)
      cat(paste("\n", dfName, " has been written to the global environment", sep=""))
      
      filenamePara <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_df_para.csv", sep=""))
      write.csv(paramPara, file=filenamePara, row.names=FALSE)
      
      cat(paste("\nSaving file: ", filenamePara,  sep=""))
      cat(paste("\n", projectName, "_para_df.csv can be opened in MS Excel.\n",  sep=""))
    }
  }
}

.findFoG <- function(i, data, ML, ML2, stand, clearHaloStand, dotedge = 3.4, maxDist = 35, standardLoc = 2.5, standType = standType){
  
  # this helper function is mostly taken from the CRAN legacy version
  # changes from that version are marked with comments
  
  startX <- which(data[[i]][,1] > dotedge+0.5)[1]
  stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  
  # CHANGE: the CRAN legacy version does not account for standType
  # this version of the code uses code from the GitHub legacy version to add this ability (also used in maxLik)
  if (standType == 'one') {
    data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand
  } else {
    data[[i]]$x <- data[[i]]$x -min(data[[i]]$x)
  }
  
  data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
  
  xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  yy<- .curveDoubleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx)
  ploty <- data[[i]]$x
  ploty[ploty < 0] <-0
  
  if (!is.null(ML)) {
    asym <- (ML[[i]]$par[1]+min(data[[i]]$x))
  } else {
    asym <- ML2[[i]]$par[1]+min(data[[i]]$x)
  }
  
  xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  yy <- (yy+min(data[[i]]$x))
  yy[yy < 0] <- 0
  
  if(max(yy) < asym*0.8){
    x80 <- xx[which.max(yy> max(yy) * 0.8)]
    x50 <- xx[which.max(yy> max(yy) * 0.5)]
    x20 <- xx[which.max(yy> max(yy) * 0.2)]
  }
  
  else{
    x80 <- xx[which.max(yy> asym * 0.8)]
    x50 <- xx[which.max(yy> asym * 0.5)]
    x20 <- xx[which.max(yy> asym * 0.2)]
  }
  
  if (x80 < x50) x80 <- xx[which.max(yy> yy[length(yy)] * 0.8)]
  
  if(exp(x80)>1) xx80 <- seq(log(data[[i]]$distance[1]), log(round(exp(x80))), length=200)
  else xx80 <- seq(log(data[[i]]$distance[1]), log(data[[i]]$distance[2]), length=200)
  
  if(exp(x50)>1) xx50 <- seq(log(data[[i]]$distance[1]), log(round(exp(x50))), length=200)
  else xx50 <- seq(log(data[[i]]$distance[1]), log(data[[i]]$distance[2]), length=200)
  
  if(exp(x20)>1) xx20 <- seq(log(data[[i]]$distance[1]), log(round(exp(x20))), length=200)
  else xx20 <- seq(log(data[[i]]$distance[1]), log(data[[i]]$distance[2]), length=200)
  
  # CHANGE: different curve names are used to match the updated maxLik
  yy <- .curveDoubleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx)
  yy80 <- .curveDoubleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx80)
  yy50<- .curveDoubleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx50)
  yy20<- .curveDoubleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx20)
  
  yy <- (yy+min(data[[i]]$x))
  yy[yy < 0] <- 0.1
  yy80 <- (yy80+min(data[[i]]$x))
  yy80[yy80 < 0] <- 0.1
  yy50 <- (yy50+min(data[[i]]$x))
  yy50[yy50 < 0] <- 0.1
  yy20 <- (yy20+min(data[[i]]$x))
  yy20[yy20 < 0] <- 0.1
  
  id <- order(xx)
  id80 <- order(xx80)
  id50 <- order(xx50)
  id20 <- order(xx20)
  
  # diff calculates the difference between pairs of consecutive elements of a numeric vector
  # for each data point, rollmean calculates the average of it and the point before it
  # so together, these two calculate the rough area of a rectangle
  # then all the rectangles are summed together, comparative to taking an integral
  maxFoG <- sum(diff(xx[id])*zoo::rollmean(yy[id], 2))
  maxFoG80 <- exp(x80)*max(yy80)
  maxFoG50 <- exp(x50)*max(yy50)
  maxFoG20 <- exp(x20)*max(yy20)
  
  FoG80 <- sum(diff(exp(xx80[id80]))*zoo::rollmean(yy80[id80], 2))
  FoG50 <- sum(diff(exp(xx50[id50]))*zoo::rollmean(yy50[id50], 2))
  FoG20 <- sum(diff(exp(xx20[id20]))*zoo::rollmean(yy20[id20], 2))
  
  # CHANGE: the CRAN legacy version returns actual FoG and maxFoG values and percentages are calculated in the main body
  # since the FoG percentages are the only values returned to the user (not the actual/max values), 
  # in this version, the FoG percentage is calculated here and is the only value returned (along with RAD values)
  param <- data.frame(x80 = round(exp(x80), digits=3), x50 = round(exp(x50), digits=3), x20 = round(exp(x20), digits=3) , FoG80 = round(FoG80/maxFoG80, digits=2), FoG50= round(FoG50/maxFoG50, digits=2), FoG20= round(FoG20/maxFoG20, digits=2))
  
  if (exp(param$x80)<1) 	param$x80 <- 0
  if (exp(param$x50)<1)	param$x50 <- 0
  if (exp(param$x20)<1)	param$x20 <- 0
  return(param)
}

.findNegLogParams <- function(i, data, ML2, stand, clearHaloStand, dotedge = 3.4, maxDist = 25) {
  
  startX <- which(data[[i]][, 1] > dotedge + 0.5)[1]
  stopX <- which(data[[i]][, 1] > maxDist - 0.5)[1]
  
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand
  data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
  
  xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200) 
  yy <- .curveNegLog(xx, ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3])
  
  yy <- yy + min(data[[i]]$x)
  yy[yy < 0] <- 0
  
  # the asym returned by ML2 doesn't always match up to the actual maximum values
  # so we just take the maximum of the data
  # would be good to add in some testing to see whether max(data) or asym creates more accurate parameters
  asym <- max(yy)
  minim <- min(yy)
  
  # finding the DRAD values - Radius of Disinhibition
  x80 <- round(exp(xx[which.min(yy > asym * 0.2)]), digits = 3)
  x50 <- round(exp(xx[which.min(yy > asym * 0.5)]), digits = 3)
  x20 <- round(exp(xx[which.min(yy > asym * 0.8)]), digits = 3)
  
  # return NA if the DRAD values don't exist
  # this is the only testing added so far to ensure that DRAD values are reasonable
  if (minim > 0.2 * asym) {
    x80 <- NA
  }
  if (minim > 0.5 * asym) {
    x50 <- NA
  }
  if (minim > 0.8 * asym) {
    x20 <- NA
  }
  
  param <- data.frame(x80 = x80, x50 = x50, x20 = x20)
  return(param)
}

.findParadoxicalParams <- function(i, data, ML2, stand, clearHaloStand, dotedge = 3.4, maxDist = 25) {
  
  startX <- which(data[[i]][, 1] > dotedge + 0.5)[1]
  stopX <- which(data[[i]][, 1] > maxDist - 0.5)[1]
  
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand
  data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
  
  xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200) 
  yy <- .curveParadox(xx, ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[4], ML2[[i]]$par[6], ML2[[i]]$par[7], ML2[[i]]$par[8])
  
  yy <- yy + min(data[[i]]$x)
  yy[yy < 0] <- 0
  
  # xmin is the minimum point on the curve
  # xmin2 is the minimum observed data point
  # both will be returned so that the user can see the comparison
  # for most photos, they should be within 1.5 mm of each other
  ymin <- min(yy)
  xmin <- xx[which(yy == ymin)]
  
  ymin2 <- min(data[[i]]$x)
  xmin2 <- data[[i]]$distance[which(ymin2 == data[[i]]$x)]
  
  xmin <- exp(median(xmin))
  xmin2 <- median(xmin2)
  
  param <- data.frame(curve.maxInhib = xmin, observed.maxInhib = xmin2)
  
  return(param)
  
}

# in the legacy version, this code was located in the main body of the function
# it's in a helper function now to reduce needless repetitions
addNames <- function(data, df, nameVector, typeVector, typePlace, typeName) {
  
  # the name/line columns will be columns 1/2, followed by the data
  # if the user wants to add a type column, this will be column 3
  
  if (is.logical(nameVector)) {
    # if nameVector is logical, the name column will have the full name of the photograph
    
    if (nameVector) {
      
      # in this case, the line column will just have the first chunk of the name of the photograph
      line <- unlist(lapply(names(data), function(x) strsplit(x, '_')[[1]][1]))
      df <- data.frame(name = names(data), line, df)
      
    } else {
      
      # in this case, the line column will just be labelled with a number
      line <- seq(1, length(data))
      df <- data.frame(name = names(data), line, df)
    }
    
  } else {
    # if nameVector is not logical, the name column will have the first chunk of the name of the photograph
    
    # in this case, the line column will have the desired name supplied by the user
    line <- nameVector
    names <- unlist(lapply(names(data), function(x) strsplit(x, '_')[[1]][1]))
    df <- data.frame(name = names, line = line, df)
    
  }
  
  if (typeVector) {
    
    type <- unlist(lapply(names(data), function(x) strsplit(x, "_")[[1]][typePlace]))
    df <- data.frame(df[, 1:2], type, df[, 3:ncol(df)])
    names(df)[3] <- typeName
    
  }
  
  return(df)
}
