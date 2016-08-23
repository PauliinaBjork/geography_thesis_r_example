#
# This script contains main parts of R for my Geography thesis of snow depth modelling.
# Feel free to use this, make enhancements, alternations and correct mistakes and save your version.
# I have calculated my explanatory variables in QGIS and SAGA, and i do the raster plotting also with QGIS. 
# I recommend using RStudio for writing R because it's giving help in parameters while you write, specially if you will prepare your data in R.
#
# Steps: 
# 1. Preparing the data for prediction. The (snow) measument and all potential explanatory variable values from each measurement site are set to a data frame (one row is one site)
# 2. GAM fitting (continuous variable)
# 3. Calculating MAE, RMSE, Willmott's D  
# 4. Model evaluation using Leave-one-out-cross-validation (LOOCV) method
# 5. Local indicators of spatial autocorrelation check with pubble plot
# 6. Make a prediction for the whole area with the chosen model to raster file




# 1. Preparing data for prediction.


# Set working directory so you'll only need to write file names 
setwd("D:/path_to_mythesis")

# Read the measured data and parameters calculated in QGIS and saved in .csv format, there can be several files
snow = read.csv("snow.csv", header = TRUE)
attach(snow) 

# If needed, you can also do calculations in this phase to get the values. Column in a data frame can be referred with $

snow$snow_means<-rowMeans(lumi[4:8], na.rm=TRUE) # Mean of columns number 4 to 8 and makes a new column of the result.

# Put all needed parameters into the same data frame, and give colums reasonable names (which will show in pictures)
my_parameters<-data.frame(cbind.data.frame(site=snow$site, snow_depth=snow$snow_means, parameter1=parameterfile$parameter1, parameter2=parameterfile$parameter2, ...)

# You can refer to parameters with the name only, no need to use the dataframe name and dollar if you attach it
attach(my_parameters)

# For preliminary data-analysis (correlations etc) using RCmdr is easiest, but you can also do it from command line
install.packages("Rcmdr")
library(Rcmdr)
citation("Rcmdr")



#2. GAM fitting (continuous variable)

# MGCV package, read the manual: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.selection.html
install.packages("mgcv")
library(mgcv)
citation("mgcv") # Good idea to give credits to the packet creators as you're using open source SW.

# Start experimenting
# How is your data distribution like? If it's binary you cannot use these examples (only prediction part is useful) Is is continous? If not you have to do something else
# Is it only positive? Does the histogram resemble more gaussian or Gamma?
# Identity and log link functions (note Gamma & identity requires no zeroes in explanatory variables). Do no use Inverse link if not necessary as it reverses the results and is different to interpret.
# Try different k values (2, 3, 4, ... )

GAM.1<-gam(snow_means ~ s(parameter1, k=3) +s(parameter2, k=) + ...   , family=Gamma(log), data=my_parameters)

# Do variable selection at your preferred style. 
# For example GCV value can be minimized, it's visible in the summary. Note there can be sevaral minimums.
summary(GAM.1) 
plot(GAM.1, pages=1)
gam.check(GAM.1) 
vis.gam(GAM.1) # 3D...



# 3. Calculate prediction error metrics: MAE, RMSE, Willmott's D etc based on training data.
# Residuals are difference between predicted and measured values.
# RMSE and MAE cam be used to compare models from same data. Willmott's D serves for comparison of models with different datasets. 

install.packages("hydroGOF")
library(hydroGOF)
citation("hydroGOF")

gam1_predictors<-data.frame(parameter2, parameter5, parameter6, ...) # selected parameters
gam1_predictions<-predict.gam(GAM.1, gam1_predictors, type="response")

mae_gam1<-mae(c(gam1_predictions), c(snow_means) )
mae_gam1

rmse_gam1<-rmse(c(gam1_predictions), c(snow_means) )
rmse_gam1

d_gam1<-d(c(gam1_predictions), c(snow_means))
d_gam1

# HydroGOF is not necessarily needed, you could just calculate these but i recommend using the ready functions from packet to avoid stupid errors... there are also other error metrics available
rmse_gam1<-sqrt(mean(residuals(GAM.1, type="response")^2)) #  if your brackets are in wrong place you will get something else.




# 4. Evaluation with residuals calculated with Leave-one-out-cross-validation (LOOCV) method.
# Use this with small sample size where evaluation set cannot be separated from training set. 
# For each measurement (fold=sample size), GAM is fitted without this certain measurement, and residual is calculated between the prediction and measured value.
# Similar error metrics can be calculated for LOOCV residuals as for traning set residuals. 
# If LOOCV error metrics are not much worse than training set error metrics, model is quite robust.

install.packages("gamclass")
library(gamclass)
citation("gamclass")


cvgam1 <- CVgam(GAM.1$formula, data=my_parameters, nfold = 63, seed = 10) #seed is set to get same result every time, if can be left out as well. nfold = sample size

mae_cv_gam1<-mae(c(cvgam1$fitted), c(snow_means) )
mae_cv_gam1

rmse_cv_gam1<-rmse(c(cvgam1$fitted), c(snow_means) )
rmse_cv_gam1



# 5. Local indicators of spatial autocorrelation (LISA) consists of correlations between neighbour sites. It's easiest to interpret from bubble plot graph, not looking at the values 
# Graph consists of residuals in coordinates, when you look at the bubble plot, pay attention if there are "pockets" of similar residuals.
# If you need global Moran's I, check for example ape package.

install.packages("sp")
library(sp)
citation("sp")

# Make a data frame (SpatialPointsDataFrame) of coordinates and model residuals

# Read your coordinates from file, note the sites are in the same order as measurements and response variables
coords = read.csv("coordinates.csv", header = TRUE) 
gam1residuals<-residuals(GAM.1, type = "response")
sp_data <- data.frame(resid = gam1residuals, x = coords$coord_E, y = coords$coord_N) # replace coord_N and coord_E with the correct titles of your coords dataframe 
coordinates(sp_data) <- c("x", "y")
bubble(sp_data, "resid", col = c("blue", "orange"), main = "GAM response residuals", xlab = "X-coordinates",  ylab = "Y-coordinates")


# Bubble plot in coordinates is quite informative, you can use it to other values as well, like the response variabme
snow_spdf <- data.frame(resid = snow_means, x = coords$coord_E, y = coords$coord_N)
coordinates(snow_spdf) <- c("x", "y")
bubble(snow_spdf, "resid", col = c("blue", "darkgray"), main = "Snow depth distribution", xlab = "X-coordinates",  ylab = "Y-coordinates")




# 6. Make a prediction (raster file) for the whole area with the chosen GAM
# Reading the raster-package manual & vignette is recommended

# Raster library for raster operations
install.packages("raster")
library(raster)
citation("raster")

# RGDAL is for reading and writing raster files
install.packages("rgdal")
library(rgdal)
citation("rgdal")

# Read the predictor (which are needed for the chosen GAM) rasters, note raster file names need to be exactly the same as predictor names (+ .tiff)


parameter1_raster <- raster("parameter1.tif")
parameter2_raster<- raster("parameter2.tif")
plot(parameter1) #testing

# Make a raster stack of all used predictors and check it looks ok 
rs <- stack(parameter1_raster, parameter2_raster)
rs 
plot(rs)

# Second raster did not have min & max values in summary, they can be set afterwards, this is not needed if everything was ok
parameter1_raster <- setMinMax(parameter1_raster) 

# Predict the whole raster with raster predictors
# Note than the predict method from raster package are different and in different order than e.g. predict from mgcv-package, don't mess with these!
prediction <- predict(object = rs, model = GAM.1, type='response') 
summary (prediction) # Just checking, minimum value and median were ok in my case, but max & mean is highly too large, my GAM was not perfect but predicton works ok.

# Prediction values are too many to be printed, but i can check the beginning & end 
head (predition) 
tail(prediction)

# Values can be rounded to zero decimals
prediction_rounded<-round(prediction)

# Now also frequency table can be printed
freq(prediction_rounded)

# Simple plotting of the prediction raster
plot(prediction)

# In this case, it's not very good as borders have NA values, and the scale of snow depth is too large
# If you want to make nice raster plots with R, study ggplot2 manual
# In my case, i just print this raster to geotiff-file and do printing in QGIS 
writeRaster(prediction_rounded, filename="prediction_rounded.tif", format="GTiff")


# Study area is smaller than the whole raster, you can take smaller area easily fom the bigger raster 
# Now i will also get rid of NA values near the borders
# Just give corner coordinates for the study with extent 
e <- extent(465000, 478000, 7757000, 7766000) 
prediction_small_area<-crop(prediction_rounded,e)
writeRaster(prediction_small_area, filename="prediction_small_area.tif", format="GTiff")


# If, for some weird reason, you'd like to to get prediction data to excel, 
# it can be done by first extracting data to data frame and making a csv file
# These take only the study area extent
prediction_df <- extract(prediction, e, df=TRUE, sp=TRUE) #df=TRUE sets the data to data frame (default is a matrix)


















