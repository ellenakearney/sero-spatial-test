#create ref grid
#1.step install the packages you need
rm(list=ls())
#INLA used to fit Bayesian models
getwd()
setwd("H://Documents/test data/")

library(raster)
library(sf)
library(tidyverse)


## STEP 1: Build reference image(s) and population raster-----------------
MMR <- st_read("00 Shapefiles/mmr_polbnda2_adm1_250k_mimu.shp")
MMR_region <- subset(MMR, ST_PCODE == "MMR002" |
                       ST_PCODE == "MMR003" |
                       ST_PCODE == "MMR007")
reference.image <- raster("01 Raw data/EVI_v6.2015.Annual.mean.1km.Data.tif") %>% 
  crop(MMR_region) %>% 
  mask(MMR_region)
plot(reference.image)
#turn them all NA
in.hh <- which(!is.na(getValues(reference.image)))  #pixel id for reference.image
reference.image[in.hh] <- 1  #so only the area we will sample has value 1
reference.coordinates <- coordinates(reference.image)[in.hh,]

par(mfrow=c(1,1))
plot(reference.image,legend=F)
bigN <- length(in.hh)   #number of pixels


save(reference.image, reference.coordinates, in.hh, MMR_region, bigN, 
     file = '02 Cleaned data/pred_input_1km.Rdata')

