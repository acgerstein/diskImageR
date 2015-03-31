library(diskImageR)

setwd("D:/workspace/diskImageR/vignettes")
IJMacro("vignette", projectDir= getwd(), pictureDir = file.path(getwd(), "pictures"))
plotRaw("vignette")
maxLik("vignette", clearHalo=1, savePDF=FALSE, ZOI="all")
saveMLParam("vignette")
createDataframe("vignette", clearHalo = 1, typeName="temp")
addType("vignette", typeName="rep")
