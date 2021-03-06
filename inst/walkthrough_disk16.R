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
IJMacro16("mTest1")

#OR To specify the appropriate directories without the pop-up interface:
IJMacro16("newProject",  "/path/to/projectDir", "/path/to/projectDir/photographs/")

###################################################################################################
#Use maximum likelihood to fit single and double logistic models to the data from each photograph.
###################################################################################################
maxLik("mTest1", maxDist = 25, standType = "indiv", needMap = TRUE, testInhib = TRUE)

###################################################################################################
#Create the dataframe.

createDataframe("mTest1", standType ="indiv", needMap = TRUE, addSIR = TRUE)

###################################################################################################
#Note that if the values are not matching up to your expected values based on your ImageJ values there are a couple of things to try.
###################################################################################################
 #The first is to use a different critical value for the radius. The default is a 20% reducation in growth, but a 50% reduction might be better.
createDataframe("newProject", standType ="indiv",  RADcrit = "50%",  needMap = TRUE, addSIR = TRUE)

#If this doesn't work then the next thing to try is to turn off 'testInhib' in maxLik and then try with both critical values
maxLik("newProject", maxDist = 25, standType = "indiv", needMap = TRUE, testInhib = FALSE)
createDataframe("newProject", standType ="indiv", needMap = TRUE, addSIR = TRUE)
createDataframe("newProject", standType ="indiv",  RADcrit = "50%",  needMap = TRUE, addSIR = TRUE)

#If none of that matches then I can go back in and try out different critical values of the radius to see what matches best.
