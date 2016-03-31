This package provides a quantitative, unbiased method to analyze photographs of disk diffusion assays for any microbial species/microbial drug combination. The method measures the radius of inhibition ("RAD", i.e., resistance) at three different cutoff values (80%, 50% and 20% growth inhibition) as well as two measures of tolerance, the fraction of growth achieved above RAD ("FoG"), and drug sensitivity ("slope", the rate of change from no growth to full growth).

## Installing

The package can be installed directly from GitHub using [devtools] (http://github.com/hadley/devtools)
```r
#install.packages("devtools"
devtools::install_github("acgerstein/diskImageR")

## Vignette

Since [diskImageR] requires third-party software (ImageJ), it was not possible to include the vignette in the package in the usual way. This readme file aims to cover some of the information that is present in the vignette, but the user is encouraged to download a pdf of the full vignette [here] (http://acgerstein.weebly.com/uploads/2/3/5/6/23564534/diskimager-helpfile.pdf). 

## Function overview (in the typical order of use)

* `IJMacro`: runs an ImageJ macro on the folder that contains the photograph to be analyzed [required]
* `readInExistingIJ`: used to read in existing ImageJ analyses [optional]
* `plotRaw`: used to plot the results of ImageJ analysis [optional]
* `maxLik`: maxmimum likelihood inference to fit logistic models to the data [required]
* `saveMLParam`: save the output of maximum likelihood analysis [optional]
* `createDataframe`: dataframe creation of all parameter estimates [required]
* `addType`: add a factor column to parameter estimate dataframe [optional]
* `aggregateData`: averages data from photographs of the same strain & type [optional]
* `calcMIC`: calculate the MIC from RAD values [optional]
* `readExistingDF`: read in an existing dataframe using a pop-up box [optional]
* `oneParamPlot`: plot a single resistance or tolerance parameter [optional]
* `twoParamplot`: plot resistance (RAD) and tolerance (FoG) at a specified cutoff value [optional]
* `threeParamplot`: plot resistance (RAD), tolerance (FoG), and sensitivity (slope) [optional]


## Required software

The first step of `diskImageR` analyzes the disk diffusion photographs in ImageJ, a free, public domain image processing program available for download (http://rsb.info.nih.gov/ij/download.html). On a Mac, if you do not already have Xcode you will need to download it from the Apple Developer tools (https://developer.apple.com/xcode/download/). You may be prompted to download and install other additional programs in the R console if any are required (depends on what is already on your computer). 

## Prepare plates and photographs
The analysis done by `diskImageR` will only be as good as the photographs taken of the disk assay plates. We use the Bencher Copymate II camera mounting system. In this setup there are two fluorescent lights on either side of the disk, oriented to minimize shadows on the plate in an otherwise dark room. We use the Canon Rebel T3i camera with an ISO 800, white balance "white fluorescent light", time 1/100s, picture stype "neutral", centre focused. Any camera of reasonably high quality should suffice, though the camera should always be set in manual rather than automatic mode, as the goal is to take photographs as consistently as possible. We also use the 2s timer to avoid potentially jostling the camera while taking images and/or having a hand/arm shadow in the picture. Plates should be photographed on a dark surface (we use black velvet) and plate labels are written on the side rather than the bottom of the plate. We have also used `diskImageR`  with the XXX scanner on XXX settings. 

Prior to analysis, images should be cropped close to the plate (as above). Because the analysis program automatically detects the disk based on size, it is important that no other similar-sized circles be present in the image (e.g., from letters in labels).

Once you have the set of photographs that you want to be analyzed together they should be placed in the same directory, with nothing else inside that directory. 
<b>Important!</b> If there are any other files within the directory the script will not run properly.

The photograph file naming scheme will be be carried throughout, thus care should be taken with naming photographs in a logical manner. The general format that we use is "strain_factor1_factor2_rep.jpg". This format will allow you to use a built-in function to average across replicate pictures from the same strain. Conversely, if you intend to do this separately (or not at all) the photographs can be named anything. 

###Run the ImageJ macro on the set of photographs
The first function of the package is `IJMacro()`. From each photograph, an ImageJ macro that is included in `diskImageR` will automatically open each photograph from the specified directory, determine where the disk is located on the plate, find the center of the disk, and draw 40mm radial lines out from the center of the disk every 5 degrees. For each line, the pixel intensity will be determined at many points along the line using the built-in `plot profile` macro from ImageJ. This data will be stored in the folder *ImageJ_out* on your computer, with one file for each photograph.

`IJMacro()` can be run in two different ways, either through a user-interface with pop-up boxes, or directly through the R console. At this point you will specify a project name, the main project directory, and the photograph directory.
The project name should ideally be fairly short (easy to type without typos!) and specific to the project. It must start with a letter, not a number or special character, but can otherwise be anything. The project name must always be specified with quotation marks around it (a surprisingly common error). Otherwise there will be a red error message (like <span style="color:red">Error in IJMacro(newProject) : object 'newProject' not found</span>).
The main project directory is the place where all files generated by the package will be saved within three directories: <i>ImageJ_out</i>, <i>parameter_files</i> and <i>figures</i>. A sub-directory will be created within each of these directories with the project name for organizational purposes, so that multiple different experiments/sets of analyses can be easily conducted from the same main project directory.
The photograph directory is the one used to store photographs from above (which has nothing except the photographs to be analyzed in it). 

<b> Important! </b> There can not be any spaces or special characters in any of the folder names that are in the path that lead to either the main project directory or the photograph directory. If there are an error box titled "Macro Error" will pop up and the script will not run (the red error message <span style="color:red">Error in tList[[i]]: subscript out of bounds</span> will also show up in the R console). 

The default assumption here and in all funtions is that the disk size is the standard 6mm. If you are using custom-sized disks you will need to specify that with the argument `diskDiam = X`, where X is the size of your disk in mm. This should also be specified in the functions `plotRaw()`, `maxLik()` and `createDataframe()`, discussed below. You will also need to change the argument `standardLoc` in `maxLik()` and `createDataframe(). standardLoc is a numberic value that indicates the location (on the disk) to use to standardize intensity across photographs. The position of standardLoc is a position that should theoretically have the same intensity in all photographs, i.e., the white of the disk. The defaul value (2.5mm) was chosen after testing of 6mm disks that contain some writing. If smaller disks are used standardLoc should be scaled appropriately. You can see where standardLoc falls in each photograph in \code{plotRaw} (the red dashed line when `plotStandardLoc = TRUE`). To suppress this standardization use `standardLoc = FALSE`.

To run the ImageJ macro through a user-interface with pop-up boxes: 
```r
IJMacro("newProject")
```

If you would prefer to avoid pop-up boxes you can directly specify the main project and photograph directory locations:
```r
IJMacro("newProject", projectDir= "/path/to/projectDir", photoDir = "/path/to/projectDir/photographs/")
```

If `IJMacro()` is unable to locate ImageJ a red error will pop-up with a message like <span style="color:red">/bin/sh: /Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub: No such file or directory</span>. The easiest solution is to move ImageJ to the default location or to specify the path to ImageJ with argument `ImageJLoc = "/path/to/ImageJ"`. 

<b> Important! </b> `IJMacro()` must run completely, without error, for everything downstream. In our experience this is the most likely step for errors to occur. Errors will be indicated in the R console should they arise, and will hopefully give you clues as to what the problem is if they are different than those described above. 

After `IJMacro()` has run successfully the output of the ImageJ analysis can be found in the *ImageJ_out* directory, though this is probably not particularly helpful unless you want to see the intensity calculations from each line. The information about the average line from each photograph can be found in the "averageLines.csv" file located in the *parameter_out* folder. This is the information that is used for all further analysis within `diskImageR` and may be useful for other purposes.

To access the output of the ImageJ analysis in a later R session use the function `readInExistingIJ()` (e.g., if you ran `IJMacro()` on a separate day then you want to conduct the downstream analyses). At this step you can also change the project name, you do not have to specify the same name that was used originally. If the name is changed new subdirectory folders will be created within the main project directory. This function will bring up a pop-up box to select the main project folder and select the directory that contains the existing ImageJ output files.

```r
readInExistingIJ("betterName") 
```

###Plot the output of ImageJ analysis
The optional function `plotRaw()` will create a PDF file of plots saved to the *figures* directory that show the average pixel intensity across all 72 lines from each photograph (i.e., the data that can be found in the "averageLines.csv" file). This function is a good check to see whether the analysis proceeded properly and in and of itself may be useful to visualize differences between different strains or experimental factors.

```{r, fig.width=6, fig.height=4}
plotRaw("newProject", showNum = TRUE, popUp = FALSE, savePDF = FALSE)
```

Many different arguments can be specified to influence the plots and the PDF that is generated, including the minimum and maximum x and y values (`xmin`, `xmax`, `ymin`, `ymax`), the number of plots in each row (`xplots`), the height and width of the PDF file (`height`, `width`), the point size (`cexPt`), and the size of the x- and y-axis font (`cexX`, `cexY`). As with all functions, you can type `?plotRaw` into the R console for all options and to see default values.

###Run the maximum likelihood analysis 
The next step is the function `maxLik()`, which uses maximum likelihood to find the logistic and double logistic equations that best describe the shape of the ImageJ output data. Our primary goal in curve fitting is to capture an underlying equation that fits the observed data. These data follow a characteristic "S-shape" curve, so the standard logistic equation is used where asym is the asymptote, od50 is the midpoint, and scal is the slope at od50 divided by asym/4. The midpoint from the single logistic is used to determine sensitivity.
$$
y = \frac{asym*exp(scal(x-od50))}{1+exp(scal(x-od50))}+N(0, \sigma)
$$

We often observed disk assays that deviated from the single logistic, either rising more linearly than expected at low cell density, or with an intermediate asymptote around the midpoint. To fascilitate fitting these curves, we also fit a double logistic, which allows greater flexibility. In practice, as the double logistic has extra parameters, it will always provide a closer fit to the underlying data, thus the results of this model are used to determine the resistance and tolerance parameters. 

$$
y = \frac{asymA*exp(scalA(x-od50A))}{1+exp(scalA(x-od50A))}+\frac{asymB*exp(scalB(x-od50B))}{1+exp(scalB(x-od50B))}+N(0, \sigma)
$$

Depending on the number of photographs to be analyzed, `maxLik()` can take a fair amount of time, upwards of an hour or more. This is due to the maximum likelihood fitting procedures, which determine the best fit parameters from multiple different starting values. The status is indicated by a series of dots (".") in the R console, with one dot per photograph. This procedure is the find.mle routine from the <a href="https://github.com/richfitz/diversitree">diversitree</a> package written by Richard Fitzjohn. If for some reason the procedure gets halted in the middle of `maxLik()` (e.g., computer is shut down) as long as R remains open it should resume where it left off when the computer is reactivated.

From these functions the plate background intensity is substracted off the intensity from all values; this should be common across all pictures taken at the same time. If you are using plates with different coloured base medium their photographs should be analzyed separately, as there will be a different background intensity from different plates. The background intensity is determined from the observed pixel intensity right beside the disk on a plate where there are no colonies in this area (e.g., the photograph on the left above). This must be specified by the user through the argument `clearHalo = X`, where X is the numbered location of the appropriate photograph. Photographs are always analyzed and organized in alphabetical order; the order can be determined by typing `names(newProject)` (no quotation marks around the project name) in the R console. In our experiments we tend to have at least one appropriate photograph with a clear halo beside the disk. A good practice, however, would be to always take a photograph of a blank plate with just the disk in the center to use for this purpose (and save it with a name like "a" so that it is always the first photograph in the list (i.e., `clearHalo = 1`). The (non)results from this photograph can be removed in the function `createDataframe()` below.

The output of `maxLik()` is a list that is saved to the R environment and a PDF file with one plot per photograph that shows the results of the model fitting (saved to the *figures* directory). Many aspects of this figure can be specified including the maximum y axis (`ymax`) the number of plots on the x axis (`xplots`), the height and width of the PDF file (`height`, `width`), the values of RAD to be plotted (one of `80`, `50`, `20`, or `all`) and FoG cutoff value to plot (one of `80`, `50`, or `20`). Once `maxLik()` has been run once (in a given R session), it does not need to be rerun to made adjustments to the PDF file; to make a new figure use the argument `needML = FALSE`. The default is to save only a single PDF file (i.e., to repeatedly overwrite the same file with different figure iterations), this can be supressed with the argument `overwrite = FALSE`. 

```{r,  fig.width=5, fig.height=4, tidy=TRUE}
maxLik("newProject", clearHalo=1, RAD="all", FoG=20, needML=TRUE, overwrite = TRUE, popUp = FALSE, savePDF = FALSE)
```

<b>[OPTIONAL] Save the maximum likelihood results</b>
If you are intersted in the nuts and bolts of the maximum likelihood parameters it is possible to save these results using the `saveMLParam()` function, which will save a CSV file into the *paramter_files* directory that contains parameter estimates for asym, od50, scal and sigma, as well as the log likelihood of the single and double logistic models.
 
```
saveMLParam("newProject")
```
 
###Create and save a dataframe of parameter estimates
The last required step is to run the function `createDataframe()` to create and save a dataframe with the drug response parameter estimates, using the best fit parameters from the logistic equations:

* <b>Resistance (<i>RAD</i>)</b>
	: asymA+asymB are added together to determine the maximum level of intensity achieved on each plate (= cell density). The level of resistance (radius of inhibition, RAD), is calculated by asking what x value (distance in mm) corresponds to the point where 80%, 50% and 20% reduction in growth occured (corresponding to *RAD80*, *RAD50*, and *RAD20*)
* <b>Tolerance (<i>FoG</i>)</b>
	: the `rollmean()` function from the `zoo` package is used to calculate the area under the curve from the disk edge to each RAD cutoff value. This achieved growth is then compared to the potential growth, i.e., the area of a rectangle with length and height equal to RAD. The calculated paramaters are thus the fraction of full growth in this region (*FoG80*, *FoG50*, *FoG20*).
* <b>Sensitivity (<i>slope</i>)</b>
	: the ten data points on either side of the midpoint from the single logistic equation (od50) are used to find the slope of the best fit linear model using the lm function in R, i.e., the slope at the midpoint.
	
If you have included a blank photograph to use for the background subtraction step in `maxLik()` this can be removed from the dataframe with the argument `removeClear = TRUE`. A CSV file is written to the *parameter_files* directory which can be opened in Microsoft Excel or any program that opens text files. The dataframe is also written and saved to the R console, should you wish to conduct further analyses in R.

```{r}
createDataframe("newProject", clearHalo = 1, typeName="Temp")
newProject.df
```
If you want to access this dataframe in a later R session you can do so with the function `readExistingDF("betterName")`. Any project name can be used here, not only the previous name. This file can also be loaded in standard ways (e.g., `new.df <- read.csv(file)`) though if you intend to use the functions below, you need to save it with a name that ends with ".df".

<b>[OPTIONAL] Add additional factor columns</b>
If your photograph names contain more than one factor that is important (i.e, if your files names look like: line_factor1_factor2...") you can add extra factors into the dataframe using the function `addType()`. 

```{r}
addType("newProject", typeName="rep")
newProject.df
```

###Aggregate replicate pictures
The function `aggregateData()` is used if you have done replicate disk assays on the same strain and want to calculate their average and variance. The variance function can be specified with basic R variance measures (e.g, standard deviation, `sd`), the standard error (`se`), or the coefficient of variantion (`CV`). 

For this example I am loading an existing dataset that I call `manyReps.df`. This dataset contains data for seven different lines, with twelve replicates per line, and a factor I'm interested in that has two levels. I then use `aggregateData()` to average among the 12 replicates and calculate their standard error. `aggregateData()` will save a CSV file into the *parameter_files* directory as well as a new dataframe to the console (`manyReps.ag`).

```{r}
manyReps.df <- read.csv(file.path(getwd(), "data", "manyReps_df.csv"))
head(manyReps.df)
```

```{r}
aggregateData("manyReps", replicate=c("line", "type"), varFunc="se")
manyReps.ag
```

###Calculate MIC
The function `calcMIC` is used to convert the RAD values calculated here into the typical MIC values you would acquire with a broth microdilution assay (or an Etest strip). This conversion can be based on a) existing built-in data from a nubmer of species/drug combinations (see below), b) a user-supplied slope and intercept of the linear or quadratic relationship between RAD and log2(MIC) for the species/drug combination of interest, c) a user supplied file containing MIC information from lines previously analyzed by diskImageR for RAD, or d) a user supplied file containing both RAD and MIC information. Note that for user-supplied data (c or d) the data should not already be transformed and the file should be a .CSV file containing columns labelled "MIC" and "RAD". If the user has supplied their own MIC data the function will first determine whether a linear or quadratic model provides a better fit. A figure that plots the standard curve will be saved in the file "RAD-MIC_standardCurve.pdf" in the *figures* directory and the calculated model parameters will be saved in the file "RAD-MIC_parameters.csv" in the *parameters_out* directory.  In all cases a column containing the MIC information is added to the dataframe and the dataframe is saved. 

Either a `diskImageR` dataframe (e.g., newProject.df) or aggregated dataframe (e.g., newProject.ag) can be used.


###Plot parameter results
Three related plotting functions are included with `diskImageR`. The function `oneParamPlot()`will plot any of the single parameters (argument `param` supports "RAD20", "RAD50", "RAD80", "FoG20", "FoG50", "FoG80", "slope", the default = "RAD20") while `twoParamPlot()` will plot RAD and FoG at specified cutoff values (`RAD` supports "RAD20", "RAD50", "RAD80"; `FoG` supports "FoG20", "FoG50", "FoG80"), and `threeParamPlot()` will plot RAD, FoG and slope. 

The required input for all three functions can be the dataframe from either `createDataframe()` (specified by argument `type="df"`, the default) or from `aggregateData()` (specified by argument `type="ag"`).  

`oneParamPlot()` will plot either a barplot (argument `barplot = TRUE`) or a dotplot (argument `barplot=FALSE`). In the two and three parameter plots the default is to plot FoG as a barplot and RAD and slope as a dotplot, though FoG can also be plotted as a dotplot with argument `barplot=FALSE`. 

Many aspects of these figure can be specified depending on type of dataframe and the number of parameters. Full details are provided in the accompanying package help files in R (e.g., `?oneParamPlot`).

```{r, fig.width=5, fig.height=4}
twoParamPlot("manyReps", type= "df", popUp = TRUE, savePDF =FALSE, xlabAngle = -45)
```

```{r, fig.width=5, fig.height=4, tidy=TRUE}
twoParamPlot("manyReps", type= "ag", popUp = TRUE, savePDF =FALSE, xlabAngle = -45, order = c(1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14), xlabels =paste(rep(manyReps.ag$line[1:7], each=2), rep(c("A", "B"), 7), sep="-"))
```

###Walkthrough of typical diskImageR use


```{r, eval=FALSE, tidy=TRUE}
#For all functions type ?functionName to bring up a help file and to see current argument default values.

###Run this function only the first time you want to run diskImageR (required to load R packages from github):

#install the devtools package. 
install.packages("devtools")

###Run the following functions everytime to use diskImageR:

#load devtools
library(devtools)

#install the diskImageR package
install_github("acgerstein/diskImageR", build_vignettes = FALSE) 

#load diskImageR
library(diskImageR)

#Run the ImageJ component, save the output. "newProject" shhould be changed to something of your choice (and then the same name used throughout); note that the quotation marks are required.
#To use a pop-up box interface:
IJMacro("newProject")
#OR To specify the appropriate directories without the pop-up interface:
IJMacro("newProject",  "/path/to/projectDir", "/path/to/projectDir/photographs/")

#Plot the result of ImageJ analysis (averaged among 72 lines drawn outward from the center of the diffusion disk). 
plotRaw("newProject")

#Use maximum likelihood to fit single and double logistic models to the data from each photograph. "clearHalo" is used to specify a picture that has a clear zone beside the disk; the intensity at this point is subtracted from all photographs and will be most accurate when photographs are taken with equal lighting without shadows. This can be a blank plate with just a disk on it (removed in the next step). RAD and FoG arguments specify values for plotting only, and do not influence analysis.
maxLik("newProject", clearHalo=1, RAD="all", FoG ="50")

#Use the logistic models from maxLik() to calculate resistance (20%, 50% and 80% reduction in growth = RAD20, RAD50, RAD80), tolerance (fraction of growth achvied above RAD relative to potential growth = FoG20, FoG50, FoG80), and sensitivity (slope at RAD50), which are saved in a CSV file. If you used a blank plate for clearHalo, remove with argument removeClear = TRUE.
createDataframe("newProject", clearHalo = 1)

#[OPTIONAL] Calculate the mean and variance for parameter estimates across replicate pictures
aggregateData("newProject")

#[OPTIONAL[ Calculate MIC from RAD values using built-in data for a limited number of spcies/drug combinations or user-supplied data.
calcMIC("newProject")

```

###Acknowledgements
* Richard Fitzjohn: contributed the maximum likelihood function `find.mle()` (from the `diversitree` package	)
* Inbal Hecht: coded portions of `calcMIC()` and contributed a patch to make `IJMacro()` more compatible with Windows
* Sincere thanks also to Adi Ulman for the original motivation, Noa Blutraich, Gal Benron, and Alexander Rosenberg for testing many versions of the code presented here, Yoav Ram for going through the code from the entire package, and Darren Abbey and particularly Judith Berman for philosophical discussions about how best to computationally capture the biological variation observed in disk assay experiments.

###Questions, comments, feedback? 
Please contact Aleeza Gerstein, <gerst035@umn.edu>

###Updated
Last updated February 2016
