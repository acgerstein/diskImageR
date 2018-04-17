#For all functions type ?functionName to bring up a help file and to see current argument default values.

#run this only the first time you run through this workflow
install.packages("devtools")

#load diskImageR through github
library(devtools)
install_github("acgerstein/diskImageR", build_vignettes = FALSE)
library(diskImageR)

###################################################################################################
#Run the ImageJ component, save the output.
#"newProject" shhould be changed to something of your choice (and then the same name used throughout); note that the quotation marks are required.
###################################################################################################

#To use a pop-up box interface:
IJMacro("newProject")

#OR To specify the appropriate directories without the pop-up interface:
IJMacro("newProject",  "/path/to/projectDir", "/path/to/projectDir/photographs/")

###################################################################################################
#Use maximum likelihood to fit single and double logistic models to the data from each photograph.
###################################################################################################
maxLik("newProject", maxDist = 25, standType = "indiv", needMap = TRUE, testInhib = TRUE)
createDataframe("newProject", standType ="indiv")

###################################################################################################
#Quantify growth requirement
###################################################################################################
inhibGrowPts("newProject")
