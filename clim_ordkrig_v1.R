
# =================================================================================================== #
#                                                                                                     
#   Indicator-based assessment of post-fire recovery dynamics using satellite NDVI time-series        
#   Torres J., Goncalves J., Marcos B. and Honrado J.                                                 
#                                                                                                     
#                                                                                                     
#   Ordinary Kriging of precipitation and air temperature records for Portuguese weather stations
#
#      This script reads a table with three initial columns holding the UID codes [1] and 
#   the XY point coordinates for weather stations [2:3] followed by a series of columns {4,...n} 
#   with annual records for air temperature (with _TMP_ in the middle of the column names) or 
#   precipitation (_PREC_) data and performs Ordinary Kriging interpolation for several different 
#   semi-variogram models, nugget values, ranges and estimation types. 
#
#   The nugget component values to test can vary by temperature or precipitation variables.
#   Model types include the Exponential, Spherical, Gaussian and Matern.
#
#      10-fold cross-validation is used to assess performance and determine the best model, i.e., the
#   combination maximizing the R2. A map is produced and saved using the best combination for each 
#   variable. 
#      A predefined raster mask is used for setting the spatial resolution, CRS and extent of the 
#   interpolation process. In this case, we used a mask coincident with MODIS data at 250m, CRS: 
#   WGS 1984 UTM-29N for north Portugal.
#
# =================================================================================================== #


library(rgdal)
library(maptools)
library(gstat)
library(sp)

setwd("C:/Users/JG/Desktop/CLIM")


## ------------------------------------------------------------------------------- ##
## DATA ----
## Change these according to your settings
## ------------------------------------------------------------------------------- ##


# Read climatic data
clim_data <- read.csv("./DATA/clim_data.csv")

# Define a raster mask used for predicting/interpolating values
mask <- readGDAL("./DATA/mask.tif")

# Create a spatial points object with the climatic data
sp_data <- SpatialPointsDataFrame(coords = clim_data[,2:3], 
                                  data = clim_data, 
                                  proj4string = CRS("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))


## ------------------------------------------------------------------------------- ##
## PARAMETERS ----
## Change these according to your needs
## ------------------------------------------------------------------------------- ##

# Variogram models to test/compare
# This can be changed to different setups...
modTypes <- c("Exp", "Sph", "Gau", "Mat")

# Nugget component of the variogram
# This can be changed to different setups...
nuggets_temp <- c(0, 0.1, 0.25, 0.5, 1.5)
nuggets_prec <- c(5000, 10000, 15000, 20000, 30000)

# Range component of the variogram
# This can be changed to different setups...
ranges <- c(10000, 30000, 50000, 80000)

# Estimation/optimization methods (don't change this...)
optEstimTypes <- c("OLS","REML")




## RUN THE ANALYSES ----

# Initialize objects
resByVar <- list()
fitVariogs <- list()

vars <- names(sp_data)[4:ncol(sp_data)]
#vars <- "ANN_PREC_2005"

for(v in vars){
  
  # Let's you change nugget values to be tested by variable type
  # For example total precipitation is identified by '_PREC_' in column names
  #
  if(grepl("_PREC_",v))
    nuggets <- nuggets_prec
  if(grepl("_TMP_",v))
    nuggets <- nuggets_temp
  
  # Number of test rounds (parameter combinations)
  nr <- length(modTypes)*length(nuggets)*length(optEstimTypes)*length(ranges)
    
  # Initialize the matrix holding CV evaluation results
  res <- data.frame(matrix(nrow=nr,ncol=7, dimnames=
                             list(1:nr,c("model","nugget","range","psill","estimType","R2","RMSE"))))
  
  # Assemble the formula for each variable
  form <- as.formula(paste(v,"~ 1",sep=""))
  
  # Create a temporary dataset without NA's
  comp.cases.ind <- !is.na(sp_data@data[,v])
  sp_data_tmp <- subset(sp_data, comp.cases.ind)
  z.obs <- sp_data_tmp@data[,v]
  
  i<-0
  for(nugget in nuggets){
    
    for(rg in ranges){
      
      for(modType in modTypes){
        
        for(optEstimType in optEstimTypes){
          
          i<-i+1
          cat("\n",v,"(",i,"/",nr,")","|| Nugget =",nugget,"Range =",rg,
              "Model =",modType,"Estimation =",optEstimType,"\n")
          
          # Make the variogram
          mod <- vgm(model =  modType,
                     psill =  var(z.obs,na.rm = TRUE),
                     range =  rg,
                     nugget = nugget)
          
          if(optEstimType == "OLS"){
            # Variogram fitting by Ordinary Least Sqaure
            fit <- fit.variogram(variogram(form, sp_data_tmp), model=mod)
          }
          
          if(optEstimType == "REML"){
            # Fit the variogram by REML model
            fit <- fit.variogram.reml(form, sp_data_tmp, model=mod)
          }
          
          # Plot the semi-variogram
          #plot(variogram(form,sp_data_tmp),fit,main="Semi-variogram")
          
          # Save fitted semi-variogram
          fitVariogs[[paste("d",nugget,sep="")]][[paste("d",rg,sep="")]][[modType]][[optEstimType]] <- fit
          
          # Kriging 10-fold Cross Validation
          krig.cv <- krige.cv(form,sp_data_tmp, model=fit, nfold = 10)
          
          # Check if the model evaluation was OK
          if(any(is.na(krig.cv$var1.pred))){
            next
          }else{
            # Calculate error measures
            krig.cv.R2 <- as.numeric(cor.test(z.obs,krig.cv$var1.pred)$estimate)^2  # R Squared
            krig.cv.RMSE <- sqrt(mean((krig.cv$residual)^2)) # Root Mean Square Error
            # Params
            res[i,1] <- modType
            res[i,2] <- nugget
            res[i,3] <- rg
            res[i,4] <- var(z.obs,na.rm = TRUE) #psill
            res[i,5] <- optEstimType
            # Results
            res[i,6] <- krig.cv.R2
            res[i,7] <- krig.cv.RMSE
          }
        }
      }
    }
  }
    
  # Write performance results
  write.csv(res, paste("Krig_",v,".csv",sep=""), na="")
  
  # Get the best model parameters
  bestCombn <- res[which.max(res[,6]),]
  bestNugget <- paste("d",bestCombn[1,"nugget"],sep="")
  bestRange <- paste("d",bestCombn[1,"range"],sep="")
  bestMod <- bestCombn[1,"model"]
  bestEstimType <- bestCombn[1,"estimType"]
  bestFit <- fitVariogs[[bestNugget]][[bestRange]][[bestMod]][[bestEstimType]]

  # Perform kriging interpolation using the reference mask
  cat("\n\n-> Performing kriging on variable:",v,"...............\n")

  map_new <- krige(form, sp_data_tmp, model = bestFit, newdata=mask)

  # Output the krigged variable
  writeGDAL(map_new,paste("Krig_",v,".tif",sep=""))
  
  # Plot the 'best' semi-variogram
  png(filename = paste("SemiVariog_",v,".png",sep=""))
  plot(variogram(form, sp_data_tmp), bestFit, main="Semi-variogram")
  dev.off()
  
  cat("Done.\n\n")
  
}

