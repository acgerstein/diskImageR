## ----setup, message=FALSE, echo=FALSE------------------------------------
library(knitr)
# This is necessary to direct knitr to find the 
# 'data', and other directories that contain
# files needed to execute this document
# thanks to http://stackoverflow.com/a/24585750/1036500
opts_knit$set(root.dir=normalizePath('../'))

## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----, echo=FALSE--------------------------------------------------------
library("diskImageR")

## ------------------------------------------------------------------------
print(getwd())
runIJManual("vignette", projectDir= getwd(), pictureDir = file.path(system.file("pictures", package="diskImageR"), ""), imageJLoc = "loc2")
# runIJManual("vignette", projectDir= getwd(), pictureDir = file.path(getwd(), "pictures", ""), imageJLoc = "loc2")

## ------------------------------------------------------------------------
plotRaw("vignette", popUp = FALSE)

## ------------------------------------------------------------------------
maxLik("vignette", clearHalo=1, ZOI="all", needML=TRUE, popUp = FALSE)

## ------------------------------------------------------------------------
createDataframe("vignette", clearHalo = 1, typeName="temp")
vignette.df

## ------------------------------------------------------------------------
addType("vignette", typeName="rep")
vignette.df

