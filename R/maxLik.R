#' Maximum Likelihood Inference

#' @description \code{maxLik}  uses maximum likelihood to find the  equations that best describe the shape of the imageJ output data to then fit parameters that describe resistance and tolerance.

#' @param projectName the short name in use for the project.
#' @param standType either `one` to specify that a single photograph will be used for standardization purposes or `indiv` to use each photograph independently. Defaults to `one`.
#' @param clearHalo numeric value that indicates which picture should be used to represent a clear halo (i.e., the clear space beside the disk).
#' @param diskDiam the diameter of the diffusion disk in mm, defaults to 6.
#' @param maxDist a numeric value indicating the maximum distance away from the disk to be considered. Defaults to 25mm.
#' @param xplots a numeric value indicating how many plots to plot in each row, does not influence maximum likelihood fitting.
#' @param ymax a numeric value indicating the maximum y value plotted in each graph, does not influence maximum likelihood fitting.
#' @param height a numeric value indicating the height of the pdf file generated, does not influence maximum likelihood fitting.
#' @param width a numeric value indicating the width of the pdf file generated, does not influence maximum likelihood fitting.
#' @param RAD a numeric value indicating the critical level of the radius of inhibition (i.e., resistance) parameter to plot, does not influence maximum likelihood fitting. Currently only \code{RAD} = "80" (80\% reduction in growth), \code{RAD} = "50" (50\% reduction in growth), \code{RAD} = "20" (20\% reduction in growth), and \code{RAD} = "all" are supported.
#' @param FoG a numeric value indicating the the critical level of the area under the curve to plot, does not influence maximum likelihood fitting. Current only \code{FoG} = "80" (80\% reduction in growth), \code{FoG} = "50" (50\% reduction in growth), and \code{FoG} = "20" (20\% reduction in growth) are supported. Note, the critical level for FoG can be different than that chosen for RAD.
#' @param needML a logical value indicating whether the maximum likelihood results already exist in the global environment or not. If \code{\link{maxLik}} has already been run in this session then needML can be set to FALSE, which allows the user to replot the results without the need to rerun the time consuming maximum likelihood models. Defaults to TRUE.
#' @param popUp a logical value indicating whether to pop up the figure after it has been created.
#' @param nameVector either a logial value indicating whether to plot the photograph names above the graph or not or a vector the same length as the number of pictures containing the desired names. Defaults to TRUE.
#' @param overwrite a logical value indicating whether to overwrite existing figures created on the same day for the same project name.defaults to TRUE.
#' @param plotFoG a logical value indicating whether to plot the FoG or not. Defaults to TRUE.
#' @param plotParam a logical value indicating whether to save plots containing, at minimum, the fitted logistic equation and specified RAD levels to plot, but may also include the FoG \code{plotFoG} = "TRUE" or the components of the logistic equation \code{plotCompon} = "TRUE". Defaults to TRUE.
#' @param savePDF a logical value indicating whether to save a PDF file or open a new quartz. Defaults to TRUE.
#' @param plotSub allows you to plot only a subset of photographs - indicate with a vector the corresponding numeric indices of the data you wish to plot. Photographs are numbered alphabetically by name, and the photograph numbers can also be found by using the showNum option in \code{\link{plotRaw}}. Defaults to NA, which will plot data from all photographs. Note this does not affect the analysis component, all data is always analyzed.
#' @param plotCompon plots the two terms of the double logistic equation. Defaults to FALSE.
#' @param standardLoc is a numberic value that indicates the location (on the disk) to use to standardize white intensity across photographs. The position of standardLoc is a position that should theoretically have the same intensity in all photographs, i.e., the white of the disk. The default value (2.5mm) was chosen after testing of 6mm disks that contain some writing. If smaller disks are used standardLoc should be scaled appropriately. You can see where standardLoc falls in each photograph in \code{plotRaw} (the red dashed line when `plotStandardLoc = TRUE`). To suppress this standardization use `standardLoc = FALSE`
#' @param typical if TRUE, a logistic curve will be calculated for each photo. If FALSE, the function will determine whether logistic, paradoxical, or confounding is the best fit and then calculate the curve.

#' @details \code{\link{maxLik}} searches for the maximum likelihood parameter using the pixel intensity information previously determined from \code{\link{IJMacro}}. Three types of inhibition growth were identified and each uses a specific equation:
#' For standard growth:
#' y = asymA*exp(scalA(x-od50A))\(1+exp(scalA(x-od50A)))+asymB*exp(scalB(x-od50B)))\(1+exp(scalB(x-od50B)))+N(0, sigma)
#'where asymA and asymB are the two asymptotes, od50A and odB are the midpoints (of the two curves), scalA and scalB are the slopes at odA and odB divided by asymA/4 and asymB/4, respectively.
#' For paradoxical growth (pixel intensity increases around the drug disk, drops, then rises again):
#' a <- exp(slope * (x - shift))
#' b <- 1 - 4 * exp(slope * (x - shift)) 
#' c <- exp(2 * slope * (x - shift))
#' d <- (1 + exp(slope * (x - shift))) ^ 4
#' y=drop * a * (b + c) / d + height
#'Where Slope adjusts slope of the middle curve,height: adjusts the height of the sides, shift is the movement(left/right), drop adjusts the variation in height between the sides and the bottom of the middle curve.
#' For confounding growth (pixel intensity increases around drug disk and then goes down):
#' y = asym*exp(-1*scal(x-od50))\(1+exp(-1*scal(x-od50)))+N(0, sigma)
#' where asym is the asymptote, od50 is the midpoint, scal is the slopes at od50 divided by asym/4.

#' @section Important:
#' If you choose to use standType = "one", the photograph specified with \code{clearHalo} is extremely important to determine tolerance, as the intensity beside the disk for the chosen photograph is subtracted for all photographs. Choosing the photograph to be used for this purpose is the only subjective aspect of this pipeline; lighting and camera settings will determine the degre to which the hue of the plate backbground changes among different photographs. Care should be taken to ensure that plate background will be as similar as possible among different plates. Photographs are numbered alphabetically by name, and can also be found using \code{\link{plotRaw}}, showNum = TRUE. In many experiments a suitable strain will already be included, however a good practice is to always take a photograph of a blank plate with just the disk in the center to use for this purpose (and save it with a name like "a" so that it is always the first photograph in the list (i.e., `clearHalo = 1`). The (non)results from this photograph can be later removed in the function `createDataframe()`. Note that only standType = "one" is supported when typical = FALSE.

#' @section Warning:
#' Depending on the number of photographs to be analyzed, `maxLik()` can take a fair amount of time, upwards of an hour or more. This is due to the maximum likelihood fitting procedures, which determine the best fit parameters from multiple different starting values. The status is indicated by a series of dots (".") in the R console, with one dot per photograph. If for some reason the procedure gets halted in the middle of \code{maxLik()} (e.g., computer is shut down) as long as R remains open it should resume where it left off when the computer is reactivated.

#' @return When typical = TRUE: two lists, ML and ML2 are saved to the global environment.
#' When typical = FALSE: one list, ML2 is saved to the global environment.
#' In both cases, a pdf file with one plot for each photograph is saved to visualize the results of curve fitting and parameters of interest.

#' @export

#' @section References:
#' Richard G. Fitzjohn (2012) Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution. 3:1084-1092.

#' @examples
#' \dontrun{
#' maxLik("myProject", clearHalo=1)
#' maxLik("myProject", clearHalo=1, xplots = 2, height = 4, width = 6, needML = FALSE)
#' }

# two parameters removed from the GitHub legacy version: needMap and testInhib
# these parameters are used in the 16-drug photos, which are not supported here

maxLik <- function(projectName, standType ="one", clearHalo, diskDiam = 6, standardLoc = 2.5, maxDist=25, ymax=200, xplots = 4, height = 8,  width = 8, FoG=20,  RAD="all", needML = TRUE, popUp = TRUE, nameVector=TRUE, overwrite = TRUE, plotParam = TRUE, plotFoG = TRUE, savePDF= TRUE, plotSub = NA, plotCompon=FALSE, typical = TRUE) {
  
  options(warn=-1)
  if(!RAD %in% c(80, 50, 20, "all")){
    stop("Current suppported RAD values = 'all', 80, 50, 20, 5")
  }
  
  # standType = 'indiv' is not yet supported when typical = FALSE
  # this method finds the minimum pixel intensity in the area near the drug disk and uses it to standardize
  # the area near the drug disk often has heavy growth for confounding and paradoxical curves so this method can't be used
  if (!typical & standType == 'indiv') {
    stop('Individual standardization with paradoxical/confounding curves is not currently supported')
  }
  
  # the plots created will be placed into the second directory created
  # if the directories already exist, new directories will not be created, but the warnings generated will not be shown
  fileFolder <- projectName
  dir.create(file.path(getwd(), "figures"), showWarnings= FALSE)
  dir.create(file.path(getwd(), "figures", fileFolder), showWarnings= FALSE)
  
  data <- eval(parse(text=projectName))
  
  if (is.logical(nameVector)){
    # if nameVector is TRUE, the photo names will be plotted in the graph
    if (nameVector){label <- names(data)}
    
    # if the names should not be plotted, we just fill label with empty strings so nothing shows up
    else {label <- rep("", length(data))}
    
  }
  else {
    # if a vector with names is supplied, those names will be used
    label <- nameVector
  }
  
  dotedge <- diskDiam/2+0.7
  
  if(standType=="one"){	
    
    if(!(hasArg(clearHalo))){
      cont <- readline(paste("Please specify photograph number with a clear halo: ", sep=""))
      clearHalo <- as.numeric(cont)
    }
    
    # the next line will find the max pixel intensity within the drug disk
    # standardLoc gives the distance from the center of the drug disk to a point near the edge
    # for each photo, find the first point that is greater than standardLoc and take the pixel intensity at that point
    # then we find the maximum of all the photos
    dotMax <- max(sapply(data, function(x) {x[which(x[,1] > standardLoc)[1], 2]}))
    
    # to standardize each photo, subtract pixel intensity at a uniform point from dotMax found above
    standard <-c( sapply(data, function(x) {dotMax-x[which(x[,1] > standardLoc)[1], 2]}))
    
  } else {
    # we make standard NULL when standType = 'indiv' so that helper functions can be combined
    # if this is the case, the standardization for each photo will be calculated inside those helper functions
    standard <- NULL
  }
  
  # optimization calculations: returns a large list containing one list for each photo
  # the first element of each photo's list contains the parameter estimates
  # if typical = FALSE, the last element of each photo's list contains the type of curve used
  
  # unlike the legacy version, we save the ML parameters using saveRDS each time they are calculated
  # this is not done every time in the legacy version
  if (needML) {
    
    # if typical, calculate both ML and ML2, same as the legacy version
    if (typical) {
      
      cat("\nStatus of single logistic ML: ")
      ML <-lapply(c(1:length(data)), .getstatsTyp, data=data, stand = standard, dotedge=dotedge, maxDist=maxDist, maxSlope=300, standType = standType)
      
      names(ML) <- names(data)
      
      # save on the local computer and to the local environment
      assign(paste(projectName, ".ML", sep=""), ML, inherits=TRUE)
      cat(paste("\n", projectName, ".ML has been written to the global environment\n", sep=""))
      
      filename.ML <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML", sep=""))
      saveRDS(ML, file=filename.ML)
      cat(paste0("\n", projectName, ".ML has been saved to ", filename.ML))
      
      cat("\nPlease note the following step may take up to an hour depending on the number of photographs being analyzed. Don't panic.\n")
      cat("\nStatus of double logistic ML: ")
      ML2 <- lapply(c(1:length(data)), .getstats2Typ, data=data, stand = standard, dotedge=dotedge, maxDist=maxDist, maxSlope=300, standType = standType)
      
      names(ML2) <- names(data)
      
      assign(paste(projectName, ".ML2", sep=""), ML2, inherits=TRUE)
      cat(paste("\n", projectName, ".ML2 has been written to the global environment\n", sep=""))
      
      filename.ML2 <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML2", sep=""))
      saveRDS(ML2, file=filename.ML2)
      cat(paste0("\n", projectName, ".ML2 has been saved to ", filename.ML2))
      
    } else {
      
      # the single logistic is not calculated for the standard photos, have to use the legacy version for that
      cat("\nPlease note the following step may take up to an hour depending on the number of photographs being analyzed. Don't panic.\n")
      cat("\nStatus of ML: ")
      ML2 <- lapply(c(1:length(data)), .getstatsAll, data=data, stand = standard, dotedge=dotedge, maxDist=maxDist, maxSlope=300)
      
      names(ML2) <- names(data)
      
      assign(paste(projectName, ".ML2", sep=""), ML2, inherits=TRUE)
      cat(paste("\n", projectName, ".ML2 has been written to the global environment\n", sep=""))
      
      filename.ML2 <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML2", sep=""))
      saveRDS(ML2, file=filename.ML2)
      cat(paste0("\n", projectName, ".ML2 has been saved to ", filename.ML2))
    }
  }
  
  if(!needML){
    # since parameters are saved every time after running ML calculations, you shouldn't have to run them more than once for any given set of photos
    
    # only the 'typical' version will have ML, but both versions will have ML2
    if (typical) {
      MLt <- paste(projectName, ".ML", sep="")
      
      # ls() returns a vector with the names of objects in the global environment
      if(MLt %in% ls()){
        ML <- eval(parse(text=MLt))
      } else {
        # if ML estimates can't be loaded from the global environment, it will check in the working directory
        filename.ML <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML", sep=""))
        ML <- readRDS(filename.ML)
        assign(paste(projectName, ".ML", sep=""), ML, inherits=TRUE)
      }
    }
    
    MLt2 <- paste(projectName, ".ML2", sep="")
    
    if(MLt2 %in% ls()) ML2 <- eval(parse(text=MLt2))
    else {
      filename.ML2 <- file.path(getwd(), "parameter_files", projectName, paste(projectName, "_ML2", sep=""))
      ML2 <- readRDS(filename.ML2)
      assign(paste(projectName, ".ML2", sep=""), ML2, inherits=TRUE)
    }
    cat(paste("\nUsing existing ML2 results ", MLt2, sep=""))
  }
  
  if (standType == 'one') {
    
    clearHaloData <- data[[clearHalo]]
    # we must add the 0.5 to dotedge here to make the data consistent with the filtering done in .plotParam helper functions
    startX <- which(clearHaloData[,1] > dotedge+0.5)[1]
    stopX <- which(clearHaloData[,1] > maxDist - 0.5)[1]
    clearHaloData <- clearHaloData[startX:stopX, 1:2]
    clearHaloData$x <- clearHaloData$x + standard[clearHalo]
    
    # recalculate distances to be relative to the disk edge rather than the disk center
    clearHaloData$distance <- clearHaloData$distance - (dotedge+0.5)
    clearHaloStand <- clearHaloData[1,2]
    
  } else {
    
    clearHaloStand <- NULL
    
  }
  
  # both .plotParam helper functions create and save a PDF file with a graph for each photo
  # each graph will contain raw data, a fitted curve, and parameters of interest
  if (plotParam) {
    
    if (typical) {
      
      .plotParamTyp(projectName, ML , ML2, stand = standard,  clearHaloStand = clearHaloStand, standardLoc = standardLoc, ymax=ymax, dotedge = dotedge, maxDist= maxDist, xplots = xplots, height = height, width=width,  FoG=FoG, RAD=RAD, overwrite = overwrite, popUp = popUp, plotFoG = plotFoG, label=label, savePDF = savePDF, plotSub = plotSub, plotCompon=plotCompon, standType = standType)
      
    } else {
      
      .plotParamAll(projectName, ML2, stand = standard,  clearHaloStand = clearHaloStand, standardLoc = standardLoc, ymax=ymax, dotedge = dotedge, maxDist= maxDist, xplots = xplots, height = height, width=width, overwrite = overwrite, popUp = popUp, label=label, savePDF = savePDF, plotSub = plotSub, plotCompon=plotCompon, RAD)
      
    }
  }
}

# used for ML estimates when typical = TRUE
# and for the second term of the paradoxical curves
.curveSingleLog <- function(asym, ic50,scal, x) {asym*exp(scal*(x-ic50))/(1+exp(scal*(x-ic50)))}

.curveDoubleLog <- function(asym, od50, scal, asymB, od50B, scalB, x) { asym*exp(scal*(x-od50))/(1+exp(scal*(x-od50)))+asymB*exp(scalB*(x-od50B))/(1+exp(scalB*(x-od50B)))}

.curveNegLog <- function(x, asym, od50, scal) { asym*exp(-1 * scal*(x-od50))/(1+exp(-1 * scal*(x-od50)))}

.curveParadox <-  function(x, slope, height, shift, drop, asym, midpoint, scal) {
  
  a <- exp(slope * (x - shift))
  b <- 1 - 4 * exp(slope * (x - shift))
  c <- exp(2 * slope * (x - shift))
  d <- (1 + exp(slope * (x - shift))) ^ 4
  
  drop * a * (b + c) / d + height + asym*exp(scal*(x-midpoint))/(1+exp(scal*(x-midpoint)))
}

# used when plotCompon = TRUE, plots the first term of the paradoxical curve
.curveParadoxCompon1 <- function(x, slope, height, shift, drop) {
  a <- exp(slope * (x - shift))
  b <- 1 - 4 * exp(slope * (x - shift))
  c <- exp(2 * slope * (x - shift))
  d <- (1 + exp(slope * (x - shift))) ^ 4
  
  drop * a * (b + c) / d + height
}

# combining .getstatsLog and .getstatsLogIndiv from the legacy version: transfer standType argument and set stand = NULL as the default to avoid error
# also taking the default value for maxSlope from .getstatsLogIndiv as it is larger 
# other than standType adjustments, this is the same as the legacy version

# .getstatsTyp is called when typical = TRUE: will calculate single logistic parameter estimates
.getstatsTyp <- function(i, data, stand = NULL, dotedge=dotedge, maxDist=maxDist, maxSlope=300, standType = standType){
  
  cat(".")
  # discrepancy here: in the legacy version, some helper functions add 0.5 to dotedge, some don't
  # for consistency, here we have opted to use +0.5 in every helper function
  startX <- which(data[[i]][,1] > dotedge+0.5)[1]
  stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  data[[i]] <- subset(data[[i]], data[[i]]$x != "NA")
  
  # adjust code to account for the standType
  if (standType == 'one') {
    data[[i]]$x <- data[[i]]$x+ stand[i] -min(data[[i]]$x+stand[i])
  } else {
    data[[i]]$x <- data[[i]]$x -min(data[[i]]$x[1:20])
  }
  
  data[[i]]$x[data[[i]]$x < 0] <- 0
  # adjust the distance to make it relative to the disk edge rather than center
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
  mlpointA<-find.mle(sumsquares.fit,par.tryA, method="subplex",upper=upper,lower=lower,control=list(maxit=50000))
  mlpointB<-find.mle(sumsquares.fit,par.tryB,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  mlpointC<-find.mle(sumsquares.fit,par.tryC,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  mlpointD<-find.mle(sumsquares.fit,par.tryD,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  
  mlpoint <- if (mlpointA$lnLik>mlpointB$lnLik) mlpointA else mlpointB
  mlpoint <- if (mlpointC$lnLik>mlpoint$lnLik) mlpointC else mlpoint
  mlpoint <- if (mlpointD$lnLik>mlpoint$lnLik) mlpointD else mlpoint
  mlpoint
}

# again, adding a parameter to account for the standardization type

# .getstats2Typ is called when typical = TRUE: will calculate double logistic parameter estimates
.getstats2Typ <- function(i, data, stand = NULL, dotedge=dotedge, maxDist=maxDist, maxSlope=300, standType = standType){
  
  cat(".")
  startX <- which(data[[i]][,1] > dotedge+0.5)[1]
  stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  data[[i]] <- subset(data[[i]], data[[i]]$x != "NA")
  
  if (standType == 'one') {
    data[[i]]$x <- data[[i]]$x+ stand[i] -min(data[[i]]$x+stand[i])
  } else {
    data[[i]]$x <- data[[i]]$x -min(data[[i]]$x[1:20])
  }
  
  data[[i]]$x[data[[i]]$x < 0] <- 0
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
  
  # since the double logistic 7 parameters (rather than 4 like the single logistic) we use more sets of initial parameter estimates
  # this creates the longer runtime from maxLik and the reason for the warning printed
  
  # adjusting scalA/B, sigma, od50A/B
  par.tryA <-c(asym = 0.9*highOD, od50 = log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2, asymB = 0.9*highOD, od50B = log(maxDist)/4, scalB = maxSlope*0.01)
  par.tryB <-c(asym = 0.9*highOD, od50 = log(maxDist)/4, scal = maxSlope*0.1, sigma =  0.2, asymB = 0.9*highOD, od50B = log(maxDist)/4, scalB = maxSlope*0.1)
  par.tryC<-c(asym = 0.9*highOD, od50 = log(maxDist)/2, scal =  maxSlope*0.01, sigma = 0.1, asymB = 0.9*highOD,od50B = log(maxDist)/2, scal =  maxSlope*0.01)
  par.tryD<-c(asym = 0.9*highOD, od50 = log(maxDist)/2, scal =  maxSlope*0.1, sigma = 0.1, asymB = 0.9*highOD,od50B = log(maxDist)/2, scalB =  maxSlope*0.1)
  # adjusting asym and scalA/B
  par.tryE <-c(asym = 0.5*highOD, od50 =  log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2, asymB = 0.7*highOD, od50B =  log(maxDist)/2, scalB = maxSlope*0.01)
  par.tryF <-c(asym = 0.5*highOD, od50 =  log(maxDist)/4, scal = maxSlope*0.1, sigma =  0.2, asymB = 0.7*highOD, od50B =  log(maxDist)/2, scalB = maxSlope*0.01)
  par.tryG <-c(asym = 0.5*highOD, od50 =  log(maxDist)/4, scal = maxSlope*0.1, sigma =  0.2, asymB = 0.7*highOD, od50B =  log(maxDist)/2, scalB = maxSlope*0.1)
  par.tryH <-c(asym = 0.5*highOD, od50 =  log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2, asymB = 0.7*highOD, od50B =  log(maxDist)/2, scalB = maxSlope*0.01)
  
  mlpoint<-c()
  mlpointA<-find.mle(sumsquares.fit,par.tryA, method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  mlpointB<-find.mle(sumsquares.fit,par.tryB,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  mlpointC<-find.mle(sumsquares.fit,par.tryC,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  mlpointD<-find.mle(sumsquares.fit,par.tryD,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  mlpointE<-find.mle(sumsquares.fit,par.tryE,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  mlpointF<-find.mle(sumsquares.fit,par.tryF,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  mlpointG<-find.mle(sumsquares.fit,par.tryG,method="subplex",upper=upper,lower=lower,control=list(maxit=50000))
  mlpointH<-find.mle(sumsquares.fit,par.tryH,method="subplex",upper=upper,lower=lower, control=list(maxit=50000))
  
  mlpoint <- if (mlpointA$lnLik>mlpointB$lnLik) mlpointA else mlpointB
  mlpoint <- if (mlpointC$lnLik>mlpoint$lnLik) mlpointC else mlpoint
  mlpoint <- if (mlpointD$lnLik>mlpoint$lnLik) mlpointD else mlpoint
  mlpoint <- if (mlpointE$lnLik>mlpoint$lnLik) mlpointE else mlpoint
  mlpoint <- if (mlpointF$lnLik>mlpoint$lnLik) mlpointF else mlpoint
  mlpoint <- if (mlpointG$lnLik>mlpoint$lnLik) mlpointG else mlpoint
  mlpoint <- if (mlpointH$lnLik>mlpoint$lnLik) mlpointH else mlpoint
  mlpoint
}

# .getstatsAll is called when typical = FALSE: will calculate and choose between the three types of curves
# only standType = 'one' is supported here
# there is only one helper function for typical = FALSE that does optimization
# there are no single logistic parameter estimates for logistic curves, it only returns double logistic parameter estimates
.getstatsAll <- function(i, data, stand, dotedge=dotedge, maxDist=maxDist, maxSlope=300){
  cat(".")
  
  startX <- which(data[[i]][,1] > dotedge+0.5)[1]
  stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  data[[i]] <- subset(data[[i]], data[[i]]$x != "NA")
  
  data[[i]]$x <- data[[i]]$x+ stand[i] -min(data[[i]]$x+stand[i])
  data[[i]]$x[data[[i]]$x < 0] <- 0
  
  data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
  data[[i]]$distance <- log(data[[i]]$distance)
  
  # double logistic curve
  sumsquares.fit1 <- function(theta){
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
  
  # negative logistic curve
  sumsquares.fit2 <- function(theta){
    asym<-theta[[1]]
    od50<-theta[[2]]
    scal<-theta[[3]]
    sigma<-theta[[4]]
    y<-data[[i]]$x
    x<-data[[i]]$distance
    res <- dnorm(y, (asym*exp(-1 * scal*(x-od50))/(1+exp(-1 * scal*(x-od50)))), sigma, log= T)
    sum(res)
  }
  
  # paradoxical curve
  sumsquares.fit3 <- function(theta){
    slope<-theta[[1]]
    height<-theta[[2]]
    shift<-theta[[3]]
    drop<-theta[[4]]
    sigma<-theta[[5]]
    asym<-theta[[6]]
    midpoint<-theta[[7]]
    scal <- theta[[8]]
    y<-data[[i]]$x
    x<-data[[i]]$distance
    
    a <- exp(slope * (x - shift))
    b <- 1 - 4 * exp(slope * (x - shift))
    c <- exp(2 * slope * (x - shift))
    d <- (1 + exp(slope * (x - shift))) ^ 4
    e <- asym*exp(scal*(x-midpoint))/(1+exp(scal*(x-midpoint)))
    
    res <- dnorm(y, drop * a * (b + c) / d + height + e, sigma, log= T)
    sum(res)
  }
  
  lowOD <- min(data[[i]]$x)
  highOD <- quantile(data[[i]]$x, 0.99, na.rm = TRUE)
  lowerLog <- c(0.8*highOD, 0, 0,0, 0, 0, 0)
  upperLog <- c(highOD, maxDist, maxSlope, 10, highOD,  maxDist, maxSlope)
  lowerNegLog <- c(0.8*highOD, 0, 0, 0)
  upperNegLog <- c(highOD, maxDist, maxSlope, 10)
  lowerParadox <- c(0, 0.8*highOD, 0, 0, 0, 0, 0, 0)
  upperParadox <- c(maxSlope, highOD, maxDist, 1000, 10, highOD, maxDist, maxSlope)
  
  # unfortunately, due to the slow speed of the optimization calculations, 
  # we only try one set of initial parameters for each curve
  # rather than 4 or 8 like .getstatsTyp/.getstats2Typ
  par.tryLog <-c(asym = 0.9*highOD, od50 = log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2, asymB = 0.9*highOD, od50B = log(maxDist)/4, scalB = maxSlope*0.01)
  par.tryNegLog <-c(asym = 0.9*highOD, ic50 = log(maxDist)/4, scal = maxSlope*0.01, sigma =  0.2)
  par.tryParadox <-c(slope = 0.1, height = 0.8*highOD, shift = 2, drop = 200, sigma =  0.2, asym = highOD * 0.9, od50 = log(maxDist), scal = maxSlope * 0.1)
  
  mlpoint<-c()
  mlpointLog<-find.mle(sumsquares.fit1,par.tryLog, method="subplex",upper=upperLog,lower=lowerLog, control=list(maxit=50000))
  mlpointNegLog<-find.mle(sumsquares.fit2,par.tryNegLog, method="subplex",upper=upperNegLog,lower=lowerNegLog, control=list(maxit=50000))
  mlpointParadox<-find.mle(sumsquares.fit3,par.tryParadox, method="subplex",upper=upperParadox,lower=lowerParadox, control=list(maxit=50000))
  
  # occasionally, curves will have higher lnLik for the paradoxical curves as it has more parameters
  # adding 20 to the Log and NegLog curves generally fixes this problem and ensures the right curve is chosen
  mlpoint <- if (mlpointLog$lnLik > mlpointNegLog$lnLik) mlpointLog else mlpointNegLog
  mlpoint <- if (mlpointParadox$lnLik > mlpoint$lnLik + 20) mlpointParadox else mlpoint
  
  df1 <- c(log = mlpointLog$lnLik + 20, negLog = mlpointNegLog$lnLik + 20, para = mlpointParadox$lnLik)
  
  # make sure the type of curve is returned so that it can be plotted properly
  df2 <- c('Log', 'NegLog', 'Para')
  maxim <- df2[which.max(df1)]
  mlpoint <- append(mlpoint, maxim)
  names(mlpoint[[length(mlpoint)]]) <- 'type'

  mlpoint
}

# .plotParamTyp is called when typical = TRUE: will put together a PDF of graphs with raw data, fitted curves, and parameters of interest
.plotParamTyp <- function(projectName, ML , ML2, stand = NULL,  clearHaloStand = NULL, standardLoc = 2.5, ymax=200, dotedge = 3.4, maxDist= 40, xplots = 4, height = 10, width=7,  FoG=50, RAD=50, overwrite = TRUE, popUp = TRUE, plotFoG = TRUE, label=label, savePDF = TRUE, plotSub = plotSub, plotCompon=plotCompon, standType = standType){
  
  data <- eval(parse(text=projectName))
  
  if(is.na(plotSub[1])){
    plotSub <- 1:length(data)
  }
  
  fileFolder <- projectName
  dir.create(file.path(getwd(), "figures"), showWarnings= FALSE)
  dir.create(file.path(getwd(), "figures", fileFolder), showWarnings= FALSE)
  
  # naming of files here is consistent with the legacy version
  t <- file.path("figures", projectName , paste(projectName, "_FoG.pdf", sep=""))
  if (!overwrite){
    if (file.exists(t)){
      t <- file.path("figures", projectName , paste(projectName, "_FoG_2_FoG", FoG, "_RAD", RAD, ".pdf", sep=""))
      if (file.exists(t)){
        k <- 2
        while(file.exists(t)){
          k <- k+1
          t <- file.path("figures", projectName, paste(projectName, "_FoG_", k, "_FoG", FoG, "_RAD", RAD, ".pdf", sep=""))
        }
      }
    }
  }
  
  # xplots is the number of plots in a row on one page
  # xplots and yplots will change if there aren't enough plots to fill the desired amount
  if(xplots > length(plotSub)){
    xplots <- length(plotSub)
  }
  
  # yplots is the number of plots in a column on one page
  if (ceiling(length(plotSub)/xplots) < 6) {
    yplots <- ceiling(length(plotSub)/xplots)}
  else {yplots<- 6}
  
  numpages <- ceiling(length(plotSub)/(xplots*yplots))
  
  if(savePDF){
    pdf(t, width=width, height=height)
  }

  # this just adjusts the margins of the pdf
  par(mfrow=c(yplots , xplots), mar=c(1,1,1,1), oma=c(4,5,1,1))
  
  for (k in plotSub){
    
    # this helper function adds a plot to the pdf
    .singlePlotTyp(data = data, ML = ML, ML2 = ML2, dotedge = dotedge, maxDist = maxDist, ymax = ymax, stand = stand, i = k,FoG=FoG, RAD = RAD, clearHaloStand = clearHaloStand, label=label[k], plotFoG = plotFoG, plotCompon=plotCompon, standType = standType)
    
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

# .singlePlotTyp is called when typical = TRUE: will create a plot with raw data, fitted curve, and parameters of interest for a single photo and add to an open pdf file

# removed showIC argument from the legacy version, it is not referenced again
.singlePlotTyp <- function(data, ML, ML2, stand = NULL, clearHaloStand = NULL, dotedge = 3.4, maxDist = 40, ymax = 200, FoG=50, RAD=50, i, label, plotFoG = TRUE, plotCompon=FALSE, standType = standType){
  
  startX <- which(data[[i]][,1] > dotedge+0.5)[1]
  stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  
  if (standType == 'one') {
    data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand
    asym <- (ML[[i]]$par[1]+min(data[[i]]$x))
  } else {
    data[[i]]$x <- data[[i]]$x -min(data[[i]]$x)
    asym <- ML[[i]]$par[1]
  }
  
  data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
  
  xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  yy <- .curveDoubleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx)
  
  ploty <- data[[i]]$x
  ploty[ploty < 0] <-0
  
  slope <- ML[[i]]$par[3]
  ic50 <- ML[[i]]$par[2]

  # plot the raw data points
  # the differences in this if else block are carried over from the GitHub legacy version
  if (standType == 'one') {
    
    plot(data[[i]]$distance, ploty, cex=0.7, col=grey(0.7), type="p", ylim=c(0, ymax), xlim=c(0, maxDist -dotedge), xaxt="n", yaxt="n", xlab="", ylab="")
    yyplot <- (yy+min(data[[i]]$x))
    
  } else {
    
    plot(data[[i]]$distance, ploty, cex=0.7, col=grey(0.7), type="p", ylim=c(0, ymax), xlim=c(0, maxDist), xaxt="n", yaxt="n", xlab="", ylab="")
    yyplot <- yy
    
  }
  
  axis(2, labels=FALSE)
  yyplot[yyplot < 0] <- 0
  # plot the fitted line
  points(exp(xx), yyplot, type="l", col="black", lwd=3)
  # removed two lines from .singlePlot legacy version which seem to place horizontal lines at the asymptotes
  
  useAsym <- "TRUE"
  
  # recall only certain RAD levels are supported
  yy80halo <- yyplot[which.max(yyplot> asym * 0.2)]
  yy50halo <- yyplot[which.max(yyplot> asym * 0.5)]
  yy20halo <- yyplot[which.max(yyplot> asym * 0.8)]
  
  if(yy20halo < yy50halo){
    yy20halo <- yyplot[which.max(yyplot> yyplot[length(yyplot)] * 0.8)]
    useAsym <- "FALSE"
  }
  
  xx80 <- exp(xx[which.max(yyplot> asym * 0.2)])
  xx50 <- exp(xx[which.max(yyplot> asym * 0.5)])
  xx20 <- exp(xx[which.max(yyplot> asym * 0.8)])
  
  if(useAsym == "FALSE"){
    xx20 <- exp(xx[which.max(yyplot> yyplot[length(yyplot)] * 0.8)])
  }
  
  # brought in from .singlePlot for edge cases
  if(length(xx)<1){
    xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  }
  
  # the following section of FoG if statements was only found in .singleFoG, but I didn't foresee any problems using it for both

  if(FoG==80){
    xx <- exp(xx[1:which.max(exp(xx) > xx80)-1])
  }
  
  if(FoG==50){
    xx <- exp(xx[1:which.max(exp(xx) > xx50)-1])
  }
  
  if(FoG==20){
    xx <- exp(xx[1:which.max(exp(xx) > xx20)-1])
  }

  
  if(length(xx)<1){
    xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  }
  
  yy<- .curveDoubleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], log(xx))
  
  if (standType == 'one') {
    yy <- (yy+min(data[[i]]$x))
  }
  
  yy[yy < 0] <- 0
  
  if (slope >1){
    xx2 <- c(xx[1], xx, xx[length(xx)])
    yy2 <- c(0, yy, 0)
    
    # fill the area under the curve (representing fraction of possible growth)
    # only .singleFoG plots the polygon in the legacy version, but it works here for both stand types
    if(plotFoG){
      polygon(xx2, yy2, density=15, col="red")
    }
    
    points(xx, yy, type="l", col="black", lwd=2)

    if(RAD == 20){
      points(xx20, yy20halo, col="blue4", cex=2, pch=19)
    }
    
    if(RAD ==50){
      points(xx50, yy50halo, col="blue", cex=2, pch=19)
    }
    
    if(RAD ==80){
      points(xx80, yy80halo, col="deepskyblue", cex=2, pch=19)
    }

    if(RAD=="all"){
      points(xx80, yy80halo, col="blue4", cex=1.75, pch=19)
      points(xx50, yy50halo, col="blue", cex=1.75, pch=19)
      points(xx20, yy20halo, col="deepskyblue", cex=1.75, pch=19)
    }
  }
  
    # to plot each separate term of the double logistic
    if(plotCompon){
      xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
      yy2.1<- .curveSingleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3],xx)
      yy2.2<- .curveSingleLog(ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7],xx)
      
      yy1plot <- (yy2.1 +min(data[[i]]$x))
      yy2plot <- (yy2.2 +min(data[[i]]$x))
      
      yy1plot[yy1plot < 0] <-0
      yy2plot[yy2plot < 0] <-0
      points(exp(xx), yy1plot , type="l", col="orange", lwd=2, lty=2)
      points(exp(xx), yy2plot, type="l", col="orange", lwd=2, lty=2)
    }
  
  mtext(label, side=3, cex=0.6)
}

# .plotParamAll is called when typical = FALSE: will put together a PDF of graphs with raw data, fitted curves, and parameters of interest
# very similar in function to .plotParamTyp, just calls a different .singlePlot function
.plotParamAll <- function(projectName, ML2, stand = NULL,  clearHaloStand = NULL, standardLoc = 2.5, ymax=200, dotedge = 3.4, maxDist= 25, xplots = 4, height = 10, width=7, overwrite = TRUE, popUp = TRUE, label=label, savePDF = TRUE, plotSub = plotSub, plotCompon=plotCompon, RAD){
  
  data <- eval(parse(text=projectName))
  
  if(is.na(plotSub[1])){
    plotSub <- 1:length(data)
  }
  
  fileFolder <- projectName
  dir.create(file.path(getwd(), "figures"), showWarnings= FALSE)
  dir.create(file.path(getwd(), "figures", fileFolder), showWarnings= FALSE)
  
  t <- file.path("figures", projectName , paste(projectName, "_all.pdf", sep=""))
  if (!overwrite){
    if (file.exists(t)){
      t <- file.path("figures", projectName , paste(projectName, "_all2", ".pdf", sep=""))
      if (file.exists(t)){
        k <- 2
        while(file.exists(t)){
          k <- k+1
          t <- file.path("figures", projectName, paste(projectName, "_all", k, ".pdf", sep=""))
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

  par(mfrow=c(yplots , xplots), mar=c(1,1,1,1), oma=c(4,5,1,1))
  
  for (k in plotSub){
    .singlePlotAll(data = data, ML2 = ML2, dotedge = dotedge, maxDist = maxDist, ymax = ymax, stand = stand, i = k,  clearHaloStand = clearHaloStand, label=label[k], plotCompon=plotCompon, RAD)
    
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

# .singlePlotAll is called when typical = FALSE: will create a plot with raw data, fitted curve, and parameters of interest for a single photo and add to an open pdf file
.singlePlotAll <- function(data, ML2, stand, clearHaloStand, dotedge = 3.4, maxDist = 25, ymax = 200, i, label, plotCompon=FALSE, RAD){
  
  startX <- which(data[[i]][,1] > dotedge+0.5)[1]
  stopX <- which(data[[i]][,1] > maxDist - 0.5)[1]
  
  data[[i]] <- data[[i]][startX:stopX, 1:2]
  data[[i]]$x <- data[[i]]$x + stand[i] - clearHaloStand
  data[[i]]$distance <- data[[i]]$distance - (dotedge+0.5)
  
  xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  
  # plotting significant parameters: we want 80, 50, 20% of the asymptotes for Log and NegLog
  if (ML2[[i]][length(ML2[[i]])] == 'Log') {
    
    yy<- .curveDoubleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7], xx)
    asym <- ML2[[i]]$par[1]
    
  } else if (ML2[[i]][length(ML2[[i]])] == 'NegLog') {
    
    yy <- .curveNegLog(xx, ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3]) 
    asym <- ML2[[i]]$par[1]
    
  } else {
    
    yy <- .curveParadox(xx, ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[4], ML2[[i]]$par[6], ML2[[i]]$par[7], ML2[[i]]$par[8])
  }
  
  ploty <- data[[i]]$x
  ploty[ploty < 0] <-0
  
  # plotting the raw data
  plot(data[[i]]$distance, ploty, cex=0.7, col=grey(0.7), type="p", ylim=c(0, ymax), xlim=c(0, maxDist -dotedge), xaxt="n", yaxt="n", xlab="", ylab="")
  axis(2, labels=FALSE)
  
  yyplot <- yy + min(data[[i]]$x)
  
  yyplot[yyplot < 0] <- 0
  # plotting the fitted line
  points(exp(xx), yyplot, type="l", col="black", lwd=3)
  
  if(length(xx)<1){
    xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
  }
  
  # colour coding:
    # logistic curves: blue
    # negative logistic curves: red
    # paradoxical curves: green
  
  # plotting parameters of interest: using code from .singlePlotTyp
  # parameters for log and neglog curves can be found in the same way as they both look at 80, 50, and 20% reduction from the max
  # so much of their code is the same, but they are plotted using different colours
  if (ML2[[i]][length(ML2[[i]])] == 'Log' | ML2[[i]][length(ML2[[i]])] == 'NegLog') {
    
    useAsym <- "TRUE"
    
    if (ML2[[i]][length(ML2[[i]])] == 'Log') {
      yy80halo <- yyplot[which.max(yyplot> asym * 0.2)]
      yy50halo <- yyplot[which.max(yyplot> asym * 0.5)]
      yy20halo <- yyplot[which.max(yyplot> asym * 0.8)]
      
      # checking to ensure points are in the right order
      # 80% reduction in growth should happen at a point before 50% reduction in growth
      if(yy20halo < yy50halo){
        yy20halo <- yyplot[which.max(yyplot> yyplot[length(yyplot)] * 0.8)]
        useAsym <- "FALSE"
      }
      
      xx80 <- exp(xx[which.max(yyplot> asym * 0.2)])
      xx50 <- exp(xx[which.max(yyplot> asym * 0.5)])
      xx20 <- exp(xx[which.max(yyplot> asym * 0.8)])
      
      # if the points aren't in the right order, we use the final yyplot value instead of the asymptote
      if(useAsym == "FALSE"){
        xx20 <- exp(xx[which.max(yyplot> yyplot[length(yyplot)] * 0.8)])
      }
      
      col1 <- 'blue4'
      col2 <- 'blue'
      col3 <- 'deepskyblue'
      
    } else {
      
      asym <- max(yyplot)
      minim <- min(yyplot)
      
      # finding the desired pixel intensity from the curve
      yy80halo <- yyplot[which.min(yyplot> asym * 0.2)]
      yy50halo <- yyplot[which.min(yyplot> asym * 0.5)]
      yy20halo <- yyplot[which.min(yyplot> asym * 0.8)]
      
      # finding the desired distance 
      xx80 <- exp(xx[which.min(yyplot> asym * 0.2)])
      xx50 <- exp(xx[which.min(yyplot> asym * 0.5)])
      xx20 <- exp(xx[which.min(yyplot> asym * 0.8)])
      
      # checking to make sure each point exists
      # we do not check here to make sure they are in the right order, unlike above (if curve = Log)
      if (minim > 0.2 * asym) {
        xx80 <- NULL
      }
      if (minim > 0.5 * asym) {
        xx50 <- NULL
      }
      if (minim > 0.8 * asym) {
        xx20 <- NULL
      }
      
      col1 <- 'firebrick4'
      col2 <- 'firebrick3'
      col3 <- 'firebrick1'
    }
    
    # this is where we use the test above to check if the points exist
    if(RAD == 20 | RAD == 'all' & !is.null(xx20)){
      points(xx20, yy20halo, col=col1, cex=2, pch=19)
    }
    
    if(RAD ==50 | RAD == 'all' & !is.null(xx50)){
      points(xx50, yy50halo, col=col2, cex=2, pch=19)
    }
    
    if(RAD ==80 | RAD == 'all' & !is.null(xx80)){
      points(xx80, yy80halo, col=col3, cex=2, pch=19)
    }
    
  } else if (ML2[[i]][length(ML2[[i]])] == 'Para') {
    
    # xmin is the minimum point on the curve
    ymin <- min(yyplot)
    xmin <- xx[which(yyplot == ymin)]
    
    # xmin2 is the minimum observed data point
    ymin2 <- min(data[[i]]$x)
    xmin2 <- data[[i]]$distance[which(ymin2 == data[[i]]$x)]
    
    # generally, these points will be within 1.5 mm of each other but some distances are larger
    
    # sometimes, xmin or xmin2 will be vectors if more than one point meets the criteria
    # taking the median of the vector will help to minimize the distance between the calculated and observed minimum

    points(exp(median(xmin)), ymin, col = 'chartreuse', cex = 2, pch = 19)
    points(median(xmin2), ymin2, col = 'forestgreen', cex = 2, pch = 19)
    
  }

  # we can only plot the components of the logistic and paradoxical curves
  # the negative logistic curves only have one term, so nothing is plotted for those
  if(plotCompon){
    xx <- seq(log(data[[i]]$distance[1]), log(max(data[[i]][,1])), length=200)
      
    if (ML2[[i]][length(ML2[[i]])] == 'Log'){
      
      yy2.1 <- .curveSingleLog(ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3],xx)
      yy2.2 <- .curveSingleLog(ML2[[i]]$par[5], ML2[[i]]$par[6], ML2[[i]]$par[7],xx)
        
      yy1plot <- (yy2.1 +min(data[[i]]$x))
      yy2plot <- (yy2.2 +min(data[[i]]$x))
        
      yy1plot[yy1plot <0] <-0
      yy2plot[yy2plot <0] <-0
        
      points(exp(xx), yy1plot , type="l", col="orange", lwd=2, lty=2)
      points(exp(xx), yy2plot, type="l", col="orange", lwd=2, lty=2)
        
    } else if (ML2[[i]][length(ML2[[i]])] == 'Para') {
      
      yy2.1 <- .curveParadoxCompon1(xx, ML2[[i]]$par[1], ML2[[i]]$par[2], ML2[[i]]$par[3], ML2[[i]]$par[4])
      # the second term of the paradoxical curves are logistic, so we can reuse that code
      yy2.2 <- .curveSingleLog(ML2[[i]]$par[6], ML2[[i]]$par[7], ML2[[i]]$par[8],xx)

      yy1plot <- (yy2.1 +min(data[[i]]$x))
      yy2plot <- (yy2.2 +min(data[[i]]$x))
        
      yy1plot[yy1plot <0] <-0
      yy2plot[yy2plot <0] <-0
        
      points(exp(xx), yy1plot , type="l", col="orange", lwd=2, lty=2)
      points(exp(xx), yy2plot, type="l", col="orange", lwd=2, lty=2)
    }
  }
  
  mtext(label, side=3, cex=0.6)
}
