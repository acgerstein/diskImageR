#install the devtools package (only need to do once)
install.packages("devtools")
install.packages("rmarkdown")
#load devtools
library(devtools)

#install the package and vignette (need to do everytime)
install_github("acgerstein/diskImageR", build_vignettes = TRUE) 

#load diskImageR
library(diskImageR)

#For an overview
browseVignettes("diskImageR")

#Run the ImageJ analysis component, save output
#Bring each picture in the specified directory into imageJ; find the center of the disk; draw a line and find the pixel intensity, repeat every 5 degrees (72 lines); average
IJMacro("projectName")

#Plot the result of ImageJ analysis (averaged among 72 lines draft outward from the center of the diffusion disk)
plotRaw("projectName")

#Use maximum likelihood to fit a bilogistic and single logistic model to the data from each photograph
maxLik("projectName", clearHalo=1, ZOI="all")

#Use the models to calculate resistance (20%, 50% and 80% reduction in growth) and tolerance (actual growth compared to potential growth up to each resistance point), save in a .csv file 
createDataframe("projectName", clearHalo = 1)

#[OPTIONAL] Calculate the mean and error for parameter estimates across replicate pictures
aggregateData("projectName")


