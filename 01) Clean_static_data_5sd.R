###data cleaning
rm(list=ls())
library(tidyverse)
library(raster)
getwd()
setwd("H://Documents/test data/")

###Data
pr_dat <- read_csv("02 Cleaned data/test_sg6_csp_pcr_villages.csv")
fim = data.frame()
## Step 3: Read in covariates, fill in missing pixels and approximately standardise
load('02 Cleaned data/pred_input_1km.Rdata')
cov.path = "01 Raw data/"
cov.list <- list.files(cov.path, pattern=".tif") 

covariate.names <- c("accessibility_to_cities_2015_v1.0",
                     "accessibility_to_healthcare_2019",
                     "accessibility_to_healthcare_2019_walking_only",
                     "Aridity_Index_v2.Synoptic.Overall.Data.1km.Data",
                     "chirps.v2.0.2015.Annual.sum.1km.NN",
                     "distance_to_water",
                     "EVI_v6.2015.Annual.mean.1km.Data",
                     "Global_Hybrid_Pop_v2_1km_UNAdj_2015",
                     "HydroSHEDS_TWI",
                     "LST_Day_v6.2015.Annual.mean.1km.Data",
                     "LST_DiurnalDiff_v6.2015.Annual.mean.1km.Data",
                     "LST_Night_v6.2015.Annual.mean.1km.Data",
                     "MERIT_Elevation.Synoptic.Overall.Data.1km.mean",
                     "PET_v2.Synoptic.Overall.Data.1km.Data",
                     "SRTM_SlopePCT_Corrected.Synoptic.Overall.Data.1km.Data",
                     "TCB_v6.2015.Annual.mean.1km.Data",
                     "TCW_v6.2015.Annual.mean.1km.Data",
                     "tree_fraction_any_percent_2000",
                     "TSI.Martens2.Pf.2015.Annual.Mean.1km.Data",
                     "VIIRS.DNB.v2_Clean.Background.2015.Annual.Data.1km.mean")

cov_stack <- stack(paste0(cov.path, covariate.names, ".tif"))
hist(cov_stack)
transformation.types <- c("Exponential",
                          "Exponential",
                          "Exponential",
                          "Normal",
                          "Normal",
                          "Exponential",
                          "Normal",
                          "Exponential",
                          "Normal",
                          "Normal",
                          "Normal",
                          "Normal",
                          "Exponential",
                          "Normal",
                          "Exponential",
                          "Normal",
                          "Normal",
                          "Normal",
                          "Normal",
                          "Exponential"
)

Ncovariates <- length(covariate.names)
for (i in 1:Ncovariates) {cat(covariate.names[i],"\t\t",transformation.types[i],col="\n")}  #check if correct transformation for correct cov

covariates <- list()
in.country <- which(!is.na(getValues(reference.image)))
in.country.coords=coordinates(reference.image)[which(!is.na(getValues(reference.image))),]

for (k in 1:Ncovariates) {  
  cat("Processing covariate:",covariate.names[k],"...\n")
  cov.current <- raster(paste(cov.path,covariate.names[k],".tif", sep="")) %>% crop(reference.image)
  if(res(cov.current)[1] != res(reference.image)[1]){
    cov.current <- resample(cov.current, reference.image, method='bilinear')
  }
  if (!prod(dim(cov.current)==dim(reference.image))) {stop("Mismatched dimensions!\n")}
  cov.current <- raster::getValues(cov.current)[in.country]
  if (transformation.types[k]=="Normal") {
    cov.current <- (cov.current-mean(cov.current,na.rm=T))/sd(cov.current,na.rm=T)
  }
  if (transformation.types[k]=="Exponential") {
    cov.current <- cov.current+min(cov.current[cov.current>0 & cov.current!=Inf & cov.current!=-Inf & !is.nan(cov.current)],na.rm=T)+abs(min(cov.current[cov.current>0 & cov.current!=Inf & cov.current!=-Inf & !is.nan(cov.current)],na.rm=T))
    cov.current <- qnorm(pexp(cov.current,1/mean(cov.current,na.rm=T)))
  }
  nas <- which(is.na(cov.current))
  cat('NAs = ',length(nas),'\n')
  if (length(nas)>0) {
    validp <- which(!(is.na(cov.current) | cov.current==-Inf | cov.current==Inf))
    valid.coords <- in.country.coords[validp,]
    for (i in nas) {
      nearestv <- which.min((valid.coords[,1]-in.country.coords[i,1])^2+(valid.coords[,2]-in.country.coords[i,2])^2)
      cov.current[i] <- cov.current[validp[nearestv]]
    }
  }
  infs <- which(cov.current==Inf)
  cat('Infs = ',length(infs),'\n')
  if (length(infs)>0) {
    validp <- which(!(cov.current==-Inf | cov.current==Inf))
    valid.coords <- in.country.coords[validp,]
    for (i in infs) {
      nearestv <- which.min((valid.coords[,1]-in.country.coords[i,1])^2+(valid.coords[,2]-in.country.coords[i,2])^2)
      cov.current[i] <- cov.current[validp[nearestv]]
    }
  }
  neginfs <- which(cov.current==-Inf)
  cat('-Infs = ',length(infs),'\n')
  if (length(neginfs)>0) {
    validp <- which(!(cov.current==-Inf))
    valid.coords <- in.country.coords[validp,]
    for (i in neginfs) {
      nearestv <- which.min((valid.coords[,1]-in.country.coords[i,1])^2+(valid.coords[,2]-in.country.coords[i,2])^2)
      cov.current[i] <- cov.current[validp[nearestv]]
    }
  }
  hist(cov.current, main=covariate.names[k])
  cat("range: ",range(cov.current),"\n")
  Sys.sleep(2)
  covariates[[length(covariates)+1]] <- cov.current
}



covariates_static <- do.call(rbind,covariates)
covariate.names <- unlist(strsplit(covariate.names, split = '.tif'))
rownames(covariates_static)=covariate.names
#simply to put the ones who are outliers in the bounds
covariates_static[covariates_static > 5] <- 5
covariates_static[covariates_static < -5] <- -5


for(i in 1:nrow(covariates_static)){
  cc <- reference.image
  cc[in.hh] = covariates_static[i,]
  names(cc) = rownames(covariates_static)[i]
  if(i==1){ cc_stack <- stack(cc)}else{cc_stack <- stack(cc_stack, cc)}
}

covariates_static_stack <-  cc_stack

save(covariates_static, covariates_static_stack, file="02 Cleaned data/covariates_static.Rdata")
