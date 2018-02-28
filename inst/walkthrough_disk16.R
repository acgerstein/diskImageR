#For all functions type ?functionName to bring up a help file and to see current argument default values.
#load diskImageR through github
library(devtools)
install_github("acgerstein/diskImageR", build_vignettes = FALSE)
library(diskImageR)

#Run the ImageJ component, save the output. "newProject" shhould be changed to something of your choice (and then the same name used throughout); note that the quotation marks are required.

#To use a pop-up box interface:
IJMacro16("newProject")
#OR To specify the appropriate directories without the pop-up interface:
IJMacro16("newProject",  "/path/to/projectDir", "/path/to/projectDir/photographs/")

#OPTIONAL: Plot the result of ImageJ analysis (averaged among 180 lines drawn outward from the center of the diffusion disk).
plotRaw("newProject")

#Use maximum likelihood to fit single and double logistic models to the data from each photograph.
maxLik("newProject", stand = "indiv" RAD="all", FoG ="50", needMap = TRUE)

createDataframe("newProject", stand ="indiv", needMap = TRUE, needSIR = TRUE)
