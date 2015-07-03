#install the devtools package (only need to do once)
install.packages("devtools")

#load devtools
library(devtools)

#install the package and vignette (need to do everytime)
install_github("acgerstein/diskImageR", build_vignettes = TRUE) 

#load diskImageR
library(diskImageR)

#For an overview
browseVignettes("diskImageR")

#Run the ImageJ analysis component, save output
IJMacro("projectName"))

#Plot the result of ImageJ analysis (averaged among 72 lines draft outward from the center of the diffusion disk)
plotRaw("projectName")

#Use maximum likelihood to fit a bilogistic and single logistic model to the data from each photograph
maxLik("projectName", clearHalo=1, ZOI="all")

#Use the models to calculate resistance (20%, 50% and 80% reduction in growth) and tolerance (actual growth compared to potential growth up to each resistance point), save in a .csv file 
createDataframe("projectName", clearHalo = 1)

#[OPTIONAL] Calculate the mean and error for parameter estimates across replicate pictures
aggregateData("projectName")


