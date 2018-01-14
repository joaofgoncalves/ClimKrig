# ClimKrig - Ordinary Kriging of precipitation and air temperature records for Portuguese weather stations

This repository is part of the paper:
_Indicator-based assessment of post-fire recovery dynamics using satellite NDVI time-series_ (submitted to Ecological Indicators)     
By **Torres J., Gonçalves J., Marcos B. and Honrado J.**

## Description

   This script reads a table with three initial columns holding the UID codes [1] and 
   the XY point coordinates for weather stations [2:3] followed by a series of columns {4,...n} 
   with annual records for air temperature (with \_TMP\_ in the middle of the column names) or 
   precipitation (\_PREC\_) data and performs Ordinary Kriging interpolation for several different 
   semi-variogram models, nugget values, ranges and estimation types (partial-sill is kept fixed and 
   equal to variance). 

   The nugget component values to test can be set differently by temperature or precipitation variables.
   Model types include the Exponential, Spherical, Gaussian and Matern. Estimation types are Ordinary 
   Least Squares (OLS) or Restricted Maximum Likelihood (REML).

   10-fold cross-validation is used to assess performance and determine the best model, i.e., the
   combination maximizing the R<sup>2</sup>. A map is produced and saved using the best combination for each 
   variable. A predefined raster mask is used for setting the spatial resolution, CRS and extent of the 
   interpolation process. In this case, we used a mask coincident with MODIS data at 250m, CRS: 
   WGS 1984 UTM-29N for north Portugal.    

## Folder structure

- **DATA** - sample data directory
  - clim_data.csv: sample data including station IDs, XY coordinates and records for 2015 total annual precipitation and average annual temperature;
  - mask.tif: a raster layer use as a reference for interpolation.

- **RES** - interpolation results including: 
    - kriging performance (10-fold CV), 
    - 'best' variogram and 
    - interpolated variables in raster format (GeoTIFF).