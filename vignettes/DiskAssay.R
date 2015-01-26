## ----setup, message=FALSE, echo=FALSE------------------------------------
library(knitr)
# This is necessary to direct knitr to find the 
# 'data', and other directories that contain
# files needed to execute this document
# thanks to http://stackoverflow.com/a/24585750/1036500
# opts_knit$set(root.dir=normalizePath('../'))
# opts_chunk$set(fig.path = "../figures/")

## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----, echo=FALSE--------------------------------------------------------
library("diskImageR")

## ------------------------------------------------------------------------
IJMacro("vignette", projectDir= getwd(), pictureDir = file.path(getwd(), "pictures", ""))

## ----, fig.width=6, fig.height=4-----------------------------------------
plotRaw("vignette", popUp = FALSE, savePDF = FALSE)

## ----,  fig.width=5, fig.height=4----------------------------------------
maxLik("vignette", clearHalo=1, ZOI="all", needML=TRUE, popUp = FALSE, savePDF =FALSE, AUC=20)

## ------------------------------------------------------------------------
createDataframe("vignette", clearHalo = 1, typeName="temp")
vignette.df

## ------------------------------------------------------------------------
addType("vignette", typeName="rep")
vignette.df

## ------------------------------------------------------------------------
manyReps.df <- read.csv(file.path(getwd(), "data", "manyReps_df.csv"))
head(manyReps.df)

## ------------------------------------------------------------------------
aggregateData("manyReps", replicate=c("line", "type"), varFunc="se")
manyReps.ag

