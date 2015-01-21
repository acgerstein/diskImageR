---
title: "Quantitative Disk Assay"
author: "Aleeza C. Gerstein"
date: "2015-01-20"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Quantitative Disk Assay}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc} 
---

## Introduction to diskImageR

diskImageR provides a quantitative way to analyze photographs taken from disk diffusion assays, and removes the need for subjective measurement by assessing zone diameter with a ruler. This computational method measures the zone of inhibition (i.e., resistance) at three different points (80%, 50% and 20% growth inhibition) as well as two measures of tolerance, how much growth there is above the zone of inhibition (area under the curve), and the rate of change from no growth to full growth (sensitivity).


<center>
<img src="pictures/p2_30_a.JPG"  style="width: 35%; height: 35%" style="float:left," alt="" /> <img src="pictures/p1_30_a.JPG"  style="width: 35%; height: 35%" style="float:left," alt="" /> 
</center>

## Prepare plates and photographs
The analysis done by diskImageR will only be as good as your disk assay plates and the photographs you take. Plates should always be labelled on the side, not on the bottom. Care should be taken when setting up the camera to take photographs, as you want the lighting conditions to be as uniform as possible, without any shadows on the plates. Camera settings should be manual rather than automatic as much as possible. Once you have the set of photographs that you want to be analyzed together they should be placed in the same directory, with nothing else inside that directory.

Photograph naming can be used downstream to have diskImageR do a number of statistical things (e.g., averaging across replicates, caulculations of variance, t-tests). The general format is "strain_factor1_factor2_rep.pdf". Conversely, if you intend to do all the statistical analysis later, photographs can be named anything, even numbered.

Finally, photographs should be cropped carefully around the disk.

<b> Important! </b> There can not be any spaces or special characters in any of the folder names that are in the path that will lead to your pictures or the directory that you will use as the main project folder (i.e., the place where all the output files from this package will go). 

## Run the imageJ macro on the set of photographs
The first step in the diskImageR pipeline is to run the imageJ macro on the photograph directory. 

<b> Important! </b> imageJ must be installed on your computer. ImageJ is a free, public domain Java image proessing program available for download <a href="http://rsb.info.nih.gov/ij/download.html"> here</a>. Take note of the path to imageJ, as this will be needed for the first function.

From each photograph, the macro (in imageJ) will automatically determine where the disk is located on the plate, find the center of the disk, and draw 40mm lines out from the center of the disk every 5 degrees. For each line, the pixel intensity will be determined at many points along the line. This data will be stored in the folder *imageJ-out* on your computer, with one file for each photograph.

This step can be completed using one of two functions. To run the imageJ macro through a user-interface with pop-up boxes to select 
where you want the main project directory to be and where to find the location of the photograph directory:
```r
runIJ("projectName")
```

Conversely, to avoid pop-up boxes you can use the alternate function to supply the desired main project directory and photograph directory locations:


```r
runIJManual("vignette", projectDir= getwd(), pictureDir = file.path(getwd(), "pictures", ""), imageJLoc = "loc2")
```

```
## 
## Output of imageJ analyses saved in directory: /Users/acgerstein/Documents/Postdoc/Research/diskImageR/vignettes/imageJ-out/vignette/
## 
## Elements in dataframe vignette: 
## [1] "p1_30_a" "p2_30_a"
## 
```
Depending on where imageJ is located on your computer, the script may not run unless you specify the filepath. See ?runIJ for more details.

If you want to access the output of the imageJ macro in a later R session you can with
```r
readInExistingIJ("projectName") 	#can be any project name, does not have to be the same as previously used
```

### [optional] Plot the imageJ output
To plot pixel intensity from the average from all photographs use:

```r
plotRaw("vignette", savePDF=FALSE)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

## Run the maximum likelihood analysis 
The next step is to use maximum likelihood to find the logistic and double logistic equations that best describe the shape of the imageJ output data. These data follow a characteristic "S-shape" curve, so the standard logistic equation is used where asym is the asymptote, od50 is the midpoint, and scal is the slope at od50 divided by asym/4.
$$
y = \frac{asym*exp(scal(x-od50))}{1+exp(scal(x-od50))}+N(0, \sigma)
$$

We often observed disk assays that deviated from the single logistic, either rising more linearly than expected at low cell density, or with an intermediate asymptote around the midpoint. To fascilitate fitting these curves, we fit a double logistic, which allows greater flexibility. Our primary goal in curve fitting is to capture an underlying equation that fits the observed data, rather than to test what model fits better.
$$
y = \frac{asymA*exp(scalA(x-od50A))}{1+exp(scalA(x-od50A))}+\frac{asymB*exp(scalB(x-od50B))}{1+exp(scalB(x-od50B))}+N(0, \sigma)
$$

From these functions we substract off the plate background intensity from all values; this is common across all pictures taken at the same time and is determined from the observed pixel intensity on a plate with a clear halo (specified by the user). We then use the parameters identified in the logistic equations to determine the resistance parameters.

* <b>Resistance</b>
	: asymA+asymB are added together to determine the maximum level of intensity (= cell density) achieved on each plate. The level of resistance (zone of inhibition, ZOI), is calculated by asking what x value (distance in mm) corresponds to the point where 80%, 50% and 20% reduction in growth occurs (corresponding to *ZOI80*, *ZOI50*, and *ZOI20*)
* <b>Tolerance</b>
	: the 'rollmean' function from the zoo package is used to calculate the area under the curve (AUC)  in slice from the disk edge to each ZOI cutoff. This achieved growth is then compared to the potential growth, namely, the area of a rectangle with length and height equal to the ZOI. The calculated paramaters are thus the fraction of full growth in this region (*fAUC80*, *fACU50*, *fAUC20*).
* <b>Sensitivity</b>
	: the ten data points on either side of the midpoint (od50) from the single logistic equation are used to find the slope of the best fit linear model using the lm function in R.


```r
maxLik("vignette", clearHalo=1, savePDF=FALSE, ZOI="all", needML=FALSE)
```

```
## 
## Using existing ML results vignette.ML & vignette.ML2
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

### [OPTIONAL] Save the maximum likelihood results
It is possible to save the maximum likelihood results using
```
saveMLParam("vignette")
```
This will save a .csv file into the *paramter_files* directory that contains parameter estimates for asym, od50, scal and sigma, as well as the log likelihood of the single and double logistic models.
 
## Save the results 
The last required step creates and save a dataframe with the resistance parameter estimates. A .csv file is written to the *parameter_files* directory which can be opened in excel or any program that opens text files. 


```r
createDataframe("vignette", clearHalo = 1, typeName="temp")
```

```
## 
## vignette.df has been written to the global environment
## Saving file: /Users/acgerstein/Documents/Postdoc/Research/diskImageR/vignettes/parameter_files/vignette/vignette_df.csv
## vignette_df.csv can be opened in MS Excel.
```

```r
vignette.df
```

```
##      name lines temp ZOI80 ZOI50 ZOI20 fAUC80 fAUC50 fAUC20 slope
## 1 p1_30_a    p1   30    11    14    17   0.36   0.27   0.29 135.1
## 2 p2_30_a    p2   30     0    14    16   1.00   0.65   0.48 101.0
```

### [OPTIONAL] Add additional factor columns
If your photograph names contain more than one factor that is important (i.e, if your files names look like: line_factor1_factor2...") you can add extra factors into the dataframe using


```r
addType("vignette", typeName="rep")
```

```
## vignette.df has been written to the global environment
## Saving file: /Users/acgerstein/Documents/Postdoc/Research/diskImageR/vignettes/parameter_files/vignette_df.csv
## vignette_df.csv can be opened in MS Excel.
```

```r
vignette.df
```

```
##      name lines temp rep ZOI80 ZOI50 ZOI20 fAUC80 fAUC50 fAUC20 slope
## 1 p1_30_a    p1   30   a    11    14    17   0.36   0.27   0.29 135.1
## 2 p2_30_a    p2   30   a     0    14    16   1.00   0.65   0.48 101.0
```
If you want to access this dataframe in a later R session you can do so using readExistingDF("projectName"). Any project name can be used here, not only the previous name. This file can also be loaded in standard ways (e.g., temp <- read.csv(file)) though if you intend to use the functions below, you need to save it with a name that ends with ".df" (i.e., temp.df).

### [OPTIONAL] Aggregate replicate pictures
This function is useful if you have many replicate disk assays and want to calculate their average and variance. The function will calculate the standard error (se), coefficient of variantion (CV) or generic R variance measures (e.g., standard deviation, sd). 

For this example I am loading an existing dataset tha I will call "manyReps". This dataset contains data for seven different lines, with twelve replicates per line, and a factor I'm interested in that has two two levels. I can then use the function aggregateData to average among the 12 replicates and calculate their standard deviation. Note 


```r
manyReps <- read.csv("../data/manyReps_df.csv")
head(manyReps.df)
```

```
##        name line   type ZOI80 ZOI50 ZOI20 fAUC80 fAUC50 fAUC20 slope
## 1  A12_30_1  A12 levelA     0     1    10      1     NA   0.84  17.8
## 2 A12_30_10  A12 levelA     1     2     8     NA   0.39   0.65  24.9
## 3 A12_30_11  A12 levelA     0     1     9      1     NA   0.81  11.8
## 4 A12_30_12  A12 levelA     0     1    20      1     NA   0.89  14.2
## 5  A12_30_2  A12 levelA     0     1    12      1     NA   0.85  11.7
## 6  A12_30_3  A12 levelA     0     1     9      1     NA   0.76  11.9
```


```r
aggregateData("manyReps", replicate=c("line", "type"), varFunc="sd")
```

```
## 
## manyReps.ag has been written to the global environment
## 
## Saving file: /Users/acgerstein/Documents/Postdoc/Research/diskImageR/vignettes/parameter_files/manyReps/manyReps_ag.csv
## manyReps.ag can be opened in MS Excel (save as .xls file if desired)
```

```r
manyReps.ag
```

```
##    line   type ZOI80 ZOI50 ZOI20 fAUC80 fAUC50 fAUC20 slope sd.ZOI80
## 1   A12 levelA     0     1    10   0.93   0.65   0.79    18     0.29
## 2   A13 levelA    13    18    21   0.71   0.36   0.31   181     2.14
## 3   A14 levelA     9    14    19   0.56   0.39   0.38   130     1.37
## 4   A15 levelA     9    15    18   0.66   0.35   0.32   174     5.02
## 5   A16 levelA     7    14    20   0.78   0.46   0.42   117     3.60
## 6   A17 levelA    10    15    18   0.60   0.33   0.29   193     3.81
## 7   A18 levelA     3     9    12   0.78   0.51   0.42   120     2.35
## 8   A12 levelB     0     2    10   0.94   0.85   0.78    26     0.00
## 9   A13 levelB     4    16    22   0.85   0.56   0.46    74     4.77
## 10  A14 levelB     1    15    22   0.72   0.65   0.55    52     1.19
## 11  A15 levelB     6    15    20   0.74   0.44   0.42    83     5.09
## 12  A16 levelB     8    15    21   0.43   0.37   0.38    99     5.39
## 13  A17 levelB     5    13    19   0.74   0.47   0.45    85     5.98
## 14  A18 levelB     2     8    12   0.86   0.55   0.45   122     2.02
##    sd.ZOI50 sd.ZOI20 sd.fAUC80 sd.fAUC50 sd.fAUC20 sd.slope
## 1      0.62     4.23        NA        NA    0.0814     6.88
## 2      1.29     1.60    0.1725    0.0706    0.0331    22.60
## 3      0.97     1.47    0.2167    0.0706    0.0425    10.66
## 4      0.90     0.90    0.3030    0.1368    0.0744    22.90
## 5      1.09     1.82    0.2277    0.1195    0.0604    19.12
## 6      1.07     1.51    0.2756    0.1097    0.0545    17.26
## 7      0.85     1.38    0.2130    0.1065    0.0628    17.93
## 8      0.49     5.50    0.2165        NA    0.0759     8.34
## 9      1.76     2.50    0.1882    0.1839    0.1064    47.85
## 10     0.92     2.70    0.3546    0.0935    0.0623    21.44
## 11     1.56     2.71        NA    0.1106    0.0698    36.69
## 12     1.15     3.27        NA    0.1665    0.1332    39.02
## 13     5.20     1.83    0.3616        NA    0.1634    43.22
## 14     0.75     1.24    0.1568    0.0899    0.0584    18.35
```
This will also save a .csv file into the *parameter_files* directory.

##Addendum
####Forthcoming
* basic t-tests (t.test)
* basic anova (aov)
* single parameter graphics
* three parameter graphics

####Acknowledgements
Thank you to Adi Ulman for the original motivation, Noa Blutraich, Gal Benron and Alex Rosenberg for testing versions of the code presented here, and Judith Berman and Darren Abbey for philosophical discussions about how best to do this.

####Contact
Aleeza Gerstein, <gerst035@umn.edu>

