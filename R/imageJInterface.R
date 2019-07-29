#' Creates a a thin wrapper reference class for imageJ
#' Automatically verifies the locaiton of the imageJ binaries on the user's system
#'
#' @param filePath Allows the full file path to the imageJ application to be set (optional).
#' @return A reference class object
#' @author Erik White \email{erikwhite@gmail.com}
#' @export

imageJInterface <- setRefClass("imageJInterface",
  fields = list(ijFilePath = "character", memoryLimit = "numeric")
)

imageJInterface$methods(
  initialize = function(filePath = "ImageJ", memoryAllocation = 1024) {
    #Defaults to 1GB
    memoryLimit <<- memoryAllocation
    #Check the imageJ file path as soon as the object is created
    setFilePath(filePath)
  },
  setFilePath = function(filePath) {
    #Verify, store and return the location of the imageJ binaries
    #If the imageJ file path cannot be verified, return null and raise an error
    print(paste("Searching for application name or filepath: ", filePath))
    
    #Search for the imageJ application depending on operating system
    if (Sys.info()["sysname"] == "Darwin") {
      #Mac
      pathResult <- system(paste("mdfind -name ", filePath," | grep 'ImageJ.app'"), intern = TRUE)
    } else {
      #Windows or Linux
      pathResult <- Sys.which(filePath)
    }
    #Return and store the imageJ file path if found
    if (!is.na(pathResult) && nchar(pathResult) != 0) {
      ijFilePath <<- trimws(pathResult, which = "both")
      print(paste("imageJ application located at: ", ijFilePath))
      return(ijFilePath)
    } else {
      if (filePath != "ImageJ") {
        warning("The ImageJ application could not be found in the specified file path")
      } else {
        warning("The ImageJ application could not be found in the common install location on your system")
      }
      return(NULL)
    }
  },
  runScript = function(arguments) {
  #Execute a script using the bundled Java runtime
  #Returns NULL if successful, otherwise returns an error message
    #Check that the imageJ location is set
    if (!length(ijFilePath) > 0) {
      stop("The imageJ binaries have not been located. Re-initialise the imageJInterface object with the correct location for the imageJ binaries")
    }
    ijCommand <- ijFilePath
    #Running imageJ with the batch argument opens it in headless mode, runs the script and then closes it again
    arguments <- paste("-batch", arguments)
    
    #Specify correct Java runtime path
    switch(Sys.info()[['sysname']],
      Darwin = {
        #Remove .app from imageJ path
        if(compareVersion(gsub("darwin", replacement = "", R.version$os, fixed = TRUE), "10.12")) {
          #The folder structure is different on newer versions of OSX
          ijCommand <- paste0(ijCommand, "/jre/bin/java")
          arguments <- paste0("-Xmx", memoryLimit, "m -jar ", ijFilePath, "/Contents/Java/ij.jar -ijpath ", ijFilePath, " ", arguments)
        } else {
          ijCommand <- paste0(ijCommand, "/Contents/MacOS/JavaApplicationStub")
        }
      },
      Windows = {
        #It may be necessary to substitute folder path dividers
        #ijCommand <- gsub("/", "\\\\", ijCommand)
        #arguments <- gsub("/", "\\\\", arguments)
      },
      Linux  = {
        ijCommand <- paste0(ijCommand, "/jre/bin/java")
        arguments <- paste0("-Xmx", memoryLimit, "m -jar ", ijFilePath, "/ij.jar -ijpath ", ijFilePath, " ", arguments)
      }
    )
    
    tryCatch(
      expr = {
        #system2 is a more portable and flexible interface than system
        #See: https://www.rdocumentation.org/packages/base/versions/3.6.1/topics/system2
        print(paste("Command to execute:", ijCommand, arguments))
        system2(ijCommand, args = arguments, wait = TRUE)
        #message("Script executed successfully")
        return()
      },
      error = function(e){
        #message("Script failed to execute")
        return(e)
      },
      warning = function(w){
        #message("Script failed to execute")
        return(w)
      },
      finally = {}
    )
    
    #Could possibly run the imageJ script asyncronously
    #Would require the Future package: https://cran.r-project.org/web/packages/future/index.html
    #promise <- future(system2(paste(ijCommand, "-batch", arguments, sep = " ")))
    #promise %>%
    #  then(function(value) {
    #    cat("The operation completed!\n")
    #    print(value)
    #  })
  }
)