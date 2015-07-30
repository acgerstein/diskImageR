######################################################
#Do this the first time you want to run diskImageR
#(or any other package through github)
######################################################
#install the devtools package
install.packages("devtools")

######################################################
#Do this everytime to set up diskImageR
######################################################
#load devtools
library(devtools)

#install the diskImageR package
install_github("acgerstein/diskImageR", build_vignettes = FALSE) 

#load diskImageR
library(diskImageR)

#Run the ImageJ analysis component, save the output. "projectName" hshould be changed to something of your choice (and then the same name used throughout); note that the quotation marks are required.
#Bring each picture in the specified directory into imageJ; find the center of the disk; draw a line and find the pixel intensity, repeat every 5 degrees (72 lines); average
IJMacro("projectName")

#Plot the result of ImageJ analysis (averaged among 72 lines draft outward from the center of the diffusion disk). Type ?plotRaw for additional parameter options.
plotRaw("projectName")

#Use maximum likelihood to fit a bilogistic and single logistic model to the data from each photograph. "clearHalo" is used to specify a picture that has a clear halo; this is used to standardize all photographs and will be most effective when photographs are taken with equal lighting without shadows. Type ?maxLik for additional parameter options.
maxLik("projectName", clearHalo=1, ZOI="all")

#Use the models to calculate resistance (20%, 50% and 80% reduction in growth = RAD20, RAD50, RAD80), perseverence (actual growth compared to potential growth up to each resistance point = FoG20, FoG50, FoG80), and sensitivity (slope at RAD50), which are saved in a .csv file. Type ?createDataframe for additional parameter options.
createDataframe("projectName", clearHalo = 1)

#[OPTIONAL] Calculate the mean and error for parameter estimates across replicate pictures
aggregateData("projectName")


