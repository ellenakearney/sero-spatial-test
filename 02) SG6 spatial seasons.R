#SG6 fit to seasonally partioned data
#1.step install the packages you need
rm(list=ls())
getwd()
setwd("H://Documents/test data/")

list.of.packages <- c("INLA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#basic packages and parallel computing packages (add more if needed)
list.of.packages <- c("raster","malariaAtlas", "readxl","ggplot2","RColorBrewer", "ggmap", "rgdal", "rgeos","maptools", "tmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#load the packages
library(INLA)
library(raster)
library(malariaAtlas)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(sf)
library(readr)
library(mapview)
library(tmap)
library(ggplot2)


pr_dat <- read_csv("02 Cleaned data/test_sg6_villages_month.csv")
load("02 Cleaned data/covariates_static.Rdata")
load("02 Cleaned data//pred_input_1km.Rdata")

gsg6<-data.frame(id=pr_dat$id,pos_gsg6=pr_dat$sg6_pos,Nscreen = pr_dat$n, 
                 x=pr_dat$longitude, y=pr_dat$latitude,
                 month=pr_dat$month,month1=pr_dat$month+1,month_end=pr_dat$month_end,year=pr_dat$year)
gsg6<-gsg6[complete.cases(gsg6),]
head(gsg6)

gsg6 <- gsg6 %>% 
  mutate(season = dplyr::recode(month1, `1` = 1L, 
                                `2` = 2L,`3` = 2L, `4` = 2L, `5`= 2L,   `6`= 2L,  `7`= 2L,
                                `8` = 3L, `9` = 3L,  `10` = 3L, `11` = 3L,
                                `12` = 1L,`13` = 1L,
                                `14` = 2L, `15` = 2L)) 

#**********cobine the coordinate**********************
xy <- cbind(gsg6$x, gsg6$y)

dummy_coords <- data.frame(x=c(96.31,98.50, 98.126, 97.63, 97.54, 97.676, 97.69,  97.45,  98.4 , 95.876, 96.45, 96.86), 
                           y=c(17.98, 15.74, 15.39, 16.9,  18.57, 19.06 , 19.56,  19.90,  16.93, 19.217, 17.2 , 17.7))
dummy.coords <- cbind(dummy_coords$x, dummy_coords$y)
expanded_coords <- rbind(xy,dummy.coords)
#### Mesh construction
## x and y coordinates in the Response data
coords = cbind(gsg6$x, gsg6$y)
boundary = inla.nonconvex.hull(points=expanded_coords, convex=-0.1, concave=-0.2)
mesh = inla.mesh.2d(loc=expanded_coords, boundary=boundary, 
                    max.edge=c(0.1, 0.5), 
                    cutoff = 0.08)
mesh$n
plot(mesh)
points(coords, pch=20, col="red", cex=0.3)

m <- pointDistance(pr_dat[, c("longitude", "latitude")], lonlat=FALSE)
mm <- as.matrix(as.dist(m))
diag(mm) <- NA
i <- apply(mm, 1, which.min)
p <- cbind(1:nrow(mm), i) 
mm[p] 
mean(mm[p])

spde = inla.spde2.pcmatern(mesh,alpha=2 ,
                           prior.range = c(0.03,0.1), #mean distance between villages
                           prior.sigma = c(3,0.01)
)


#extract covariate at point coordinates for static covs
covariate_all<-data.frame(raster::extract(covariates_static_stack, xy))
head(covariate_all)
covs_included <- c("chirps.v2.0.2015.Annual.sum.1km.NN"    ,      
                   "distance_to_water" ,
                   "PET_v2.Synoptic.Overall.Data.1km.Data"    ,                   
                   "tree_fraction_any_percent_2000"   ,                      
                   "LST_DiurnalDiff_v6.2015.Annual.mean.1km.Data")
covariate_z <- data.frame(covariate_all[,covs_included])
head(covariate_z)
#*********cobine all data together**
gsg6cov<-cbind(gsg6,covariate_z)
head(gsg6cov)

gsg6cov_sum <- data.frame(pos_gsg6=gsg6cov$pos_gsg6,Nscreen=gsg6cov$Nscreen)
gsg6cov_mean <- data.frame(cbind(season=gsg6cov$season,x=gsg6cov$x, y=gsg6cov$y, id=gsg6cov$id,covariate_z))
#collapse down into quarters
gsg6cov <- cbind(aggregate(gsg6cov_sum,by = list(gsg6cov$season, gsg6cov$id),FUN = sum),
                 aggregate(gsg6cov_mean,by = list(gsg6cov$season, gsg6cov$id),FUN = mean)
                 
)

head(gsg6cov)
covariate_all = gsg6cov[,-1:-10]
head(covariate_all)

gsg6cov <- gsg6cov[,-5:-6]

##########################################
#run for first season - hot
gsg6cov_hot <-  filter(gsg6cov, season==1)
gsg6cov_rainy <-  filter(gsg6cov, season==2)
gsg6cov_cool <-  filter(gsg6cov, season==3)

#6. step: make predictions dataframe

covs_01 <- as.matrix(raster::extract(covariates_static_stack,reference.coordinates))
covs_02 <- as.matrix(raster::extract(covariates_static_stack,reference.coordinates))
covs_03 <- as.matrix(raster::extract(covariates_static_stack,reference.coordinates))

covs_hot <- cbind(covs_01[,covs_included])
covs_rainy <- cbind(covs_02[,covs_included])
covs_cool <- cbind(covs_03[,covs_included])

#### Mesh construction
## x and y coordinates in the Response data
sample_locations = cbind(gsg6cov_hot$x, gsg6cov_hot$y)


A <- inla.spde.make.A(mesh, sample_locations) 
indexs <- inla.spde.make.index("field", spde$n.spde)



####################################################################################


covar_z = gsg6cov_hot[,-1:-8]

stk.z = inla.stack(data=list(occurence=as.vector(gsg6cov_hot$pos_gsg6), trials=as.vector(gsg6cov_hot$Nscreen)),
                   A=list(A, 1),
                   effects=list(indexs,
                                list(
                                  # z.intercept=1,
                                  z.intercept=rep(1,nrow(gsg6cov_hot)),
                                  covariate=covar_z)),
                   tag="est.z")

formula.z = as.formula(paste("occurence ~ -1 + z.intercept + f(field, model=spde) + ",paste(names(covar_z),collapse="+")))



result_gsg6cov_hot <-   inla(formula.z, family="binomial", Ntrials=inla.stack.data(stk.z)$trials,
                             data=inla.stack.data(stk.z), verbose=FALSE,
                             control.predictor = list(compute=TRUE, A=inla.stack.A(stk.z)),
                             control.fixed=list(expand.factor.strategy='inla'),
                             control.compute = list(dic=TRUE, waic=TRUE, config=TRUE),
                             control.inla = list(int.strategy="eb"))
summary(result_gsg6cov_hot)
#############################


##

pred_covar_z <- data.frame((covs_hot))
Ap <- inla.spde.make.A(mesh = mesh, loc=reference.coordinates);dim(Ap)

mask<-reference.image; NAvalue(mask)=-9999
pred_val<-getValues(mask)
w<-is.na(pred_val)
index<-1:length(w)
index<-index[!w]
pred_locs<-xyFromCell(mask,1:ncell(mask))
pred_locs<-pred_locs[!w,]
colnames(pred_locs)<-c('longitude','latitude')
locs_pred <- pred_locs


#define number of samples 
nn = 150#(1500 for final results, 150 for fast results)

#sampling the posterior
set.seed(999)
samp = inla.posterior.sample(nn, result_gsg6cov_hot)
pred = matrix(NA,nrow=dim(Ap)[1], ncol=nn)
k = dim(pred_covar_z)[2] ## number of final covariates
for (i in 1:nn){
  field = samp[[i]]$latent[grep('field',rownames(samp[[i]]$latent)),]
  intercept = samp[[i]]$latent[grep('z.intercept',rownames(samp[[i]]$latent)),]
  beta = NULL
  for (j in 1:k){
    beta[j] = samp[[i]]$latent[grep(names(pred_covar_z)[j],rownames(samp[[i]]$latent)),]
  }
  #compute beta*covariate for each covariate
  linpred<-list()
  for (j in 1:k){
    linpred[[j]]<-beta[j]*pred_covar_z[,j]
  }
  linpred<-Reduce("+",linpred)
  lp = intercept + linpred + drop(Ap%*%field)
  ## Predicted values
  pred[,i] = plogis(lp) #for a bernoulli likelihood
  # pred[,i] = exp(lp)#for a poisson likelihood
}
pred_gsg6_mean = rowMeans(pred)
pred_gsg6_ll = apply(pred, 1, function(x) quantile(x, probs=c(0.025), na.rm=TRUE))
pred_gsg6_ul = apply(pred, 1, function(x) quantile(x, probs=c(0.975), na.rm=TRUE))
pred_gsg6_sd = apply(pred, 1, sd)
pred_gsg6_25pct = apply(pred, 1, function(x) quantile(x, probs=c(0.25), na.rm=TRUE))
pred_gsg6_75pct = apply(pred, 1, function(x) quantile(x, probs=c(0.75), na.rm=TRUE))
pred_gsg6_IQR = pred_gsg6_75pct - pred_gsg6_25pct
#Saving predictive maps as raster files
predinput=list(pred_gsg6_mean,pred_gsg6_sd,pred_gsg6_25pct,pred_gsg6_75pct,pred_gsg6_IQR,pred_gsg6_ll,pred_gsg6_ul)
prednames=as.list(c("Prob_gsg6_mean_hot_","Prob_gsg6_sd_hot_","Prob_gsg6_q25_hot_","Prob_gsg6_q75_hot_","Prob_gsg6_IQR_hot_","Prob_gsg6_ll_hot_","Prob_gsg6_ul_hot_"))
mainnames=as.list(c("PR gsg6_mean","PR gsg6_sd","PR gsg6_25th pct","PR gsg6_75th pct","gsg6_IQR","PR gsg6_ll","PR gsg6_ul"))
path_output <- paste0(getwd(),"/03 Output/gsg6/seasons")  
destfile<-paste0(path_output,'/')#make sure you have created a file OUTPUT in the bucket

out<-list()
for (j in 1:length(predinput))
{
  pred_val[!w] <- round(predinput[[j]],3)
  out[[j]] = setValues(mask, pred_val)
  writeRaster(out[[j]], paste0(destfile,prednames[[j]],',.tif'), overwrite=TRUE)
}
raster.list <- list.files(path=destfile,pattern =".tif$", full.names=TRUE)
raster.list <- raster.list[grep('hot',raster.list)]
#extract raster data
rasterls<-list()
for (i in 1:length(raster.list)){
  rasterls[[i]]<-raster(raster.list[[i]])
}
d <- brick(rasterls)
plot(d)
#add informative text before the graphs
t1 = "Predicted posterior mean probability of occurrence"
t2 = "Predicted posterior lower limit probability"
t3 = "Predicted posterior upper limit probability"
t4 = "Predicted posterior inter-quartile range"
t5 = "Predicted posterior Q25"
t6 = "Predicted posterior Q75"
t7 = "Predicted posterior standard deviation"
plot.new()
pdf.options()
pdf(file=paste0(destfile,"Prediction_probability_gsg6_states_hot.pdf"))
plot(d[[3]], main=t1)
plot(d[[2]], main=t2)
plot(d[[7]], main=t3)
plot(d[[1]], main=t4)
plot(d[[4]], main=t5)
plot(d[[5]], main=t6)
plot(d[[6]], main=t7)
dev.off();dev.off()

rb_gsg6_hot <- d



#save outputs of the analysis (beta coefficients)

fixed_table = result_gsg6cov_hot$summary.fixed#check that 
fixed_table <- fixed_table[,1:5]
fixed_table$sig <-ifelse(sign(fixed_table[,3])==sign(fixed_table[,5]),1,0)
fixed_table$or <- exp(fixed_table[,1])
fixed_table$or_lci <- exp(fixed_table[,1]+qnorm(0.025)*fixed_table[,2])
fixed_table$or_uci <- exp(fixed_table[,1]+qnorm(0.975)*fixed_table[,2])

write.csv(fixed_table, file="03 Output/gsg6/seasons/fixedeffects_gsg6_hot.csv",row.names=TRUE)


####################################################################################

########################
#### Binomial model ####
########################

#### Mesh construction
## x and y coordinates in the Response data
sample_locations = cbind(gsg6cov_rainy$x, gsg6cov_rainy$y)


A <- inla.spde.make.A(mesh, sample_locations) 
indexs <- inla.spde.make.index("field", spde$n.spde)



####################################################################################


covar_z = gsg6cov_rainy[,-1:-8]

stk.z = inla.stack(data=list(occurence=as.vector(gsg6cov_rainy$pos_gsg6), trials=as.vector(gsg6cov_rainy$Nscreen)),
                   A=list(A, 1),
                   effects=list(indexs,
                                list(
                                  # z.intercept=1,
                                  z.intercept=rep(1,nrow(gsg6cov_rainy)),
                                  covariate=covar_z)),
                   tag="est.z")

formula.z = as.formula(paste("occurence ~ -1 + z.intercept + f(field, model=spde) + ",paste(names(covar_z),collapse="+")))



result_gsg6cov_rainy <-   inla(formula.z, family="binomial", Ntrials=inla.stack.data(stk.z)$trials,
                               data=inla.stack.data(stk.z), verbose=FALSE,
                               control.predictor = list(compute=TRUE, A=inla.stack.A(stk.z)),
                               control.fixed=list(expand.factor.strategy='inla'),
                               control.compute = list(dic=TRUE, waic=TRUE, config=TRUE),
                               control.inla = list(int.strategy="eb"))

summary(result_gsg6cov_rainy)
#############################


##

pred_covar_z <- data.frame((covs_rainy))
Ap <- inla.spde.make.A(mesh = mesh, loc=reference.coordinates);dim(Ap)

mask<-reference.image; NAvalue(mask)=-9999
pred_val<-getValues(mask)
w<-is.na(pred_val)
index<-1:length(w)
index<-index[!w]
pred_locs<-xyFromCell(mask,1:ncell(mask))
pred_locs<-pred_locs[!w,]
colnames(pred_locs)<-c('longitude','latitude')
locs_pred <- pred_locs


#define number of samples 
nn = 150#(1500 for final results, 150 for fast results)

#sampling the posterior
set.seed(999)
samp = inla.posterior.sample(nn, result_gsg6cov_rainy)
pred = matrix(NA,nrow=dim(Ap)[1], ncol=nn)
k = dim(pred_covar_z)[2] ## number of final covariates
for (i in 1:nn){
  field = samp[[i]]$latent[grep('field',rownames(samp[[i]]$latent)),]
  intercept = samp[[i]]$latent[grep('z.intercept',rownames(samp[[i]]$latent)),]
  beta = NULL
  for (j in 1:k){
    beta[j] = samp[[i]]$latent[grep(names(pred_covar_z)[j],rownames(samp[[i]]$latent)),]
  }
  #compute beta*covariate for each covariate
  linpred<-list()
  for (j in 1:k){
    linpred[[j]]<-beta[j]*pred_covar_z[,j]
  }
  linpred<-Reduce("+",linpred)
  lp = intercept + linpred + drop(Ap%*%field)
  ## Predicted values
  pred[,i] = plogis(lp) #for a bernoulli likelihood
  # pred[,i] = exp(lp)#for a poisson likelihood
}
pred_gsg6_mean = rowMeans(pred)
pred_gsg6_ll = apply(pred, 1, function(x) quantile(x, probs=c(0.025), na.rm=TRUE))
pred_gsg6_ul = apply(pred, 1, function(x) quantile(x, probs=c(0.975), na.rm=TRUE))
pred_gsg6_sd = apply(pred, 1, sd)
pred_gsg6_25pct = apply(pred, 1, function(x) quantile(x, probs=c(0.25), na.rm=TRUE))
pred_gsg6_75pct = apply(pred, 1, function(x) quantile(x, probs=c(0.75), na.rm=TRUE))
pred_gsg6_IQR = pred_gsg6_75pct - pred_gsg6_25pct
#Saving predictive maps as raster files
predinput=list(pred_gsg6_mean,pred_gsg6_sd,pred_gsg6_25pct,pred_gsg6_75pct,pred_gsg6_IQR,pred_gsg6_ll,pred_gsg6_ul)
prednames=as.list(c("Prob_gsg6_mean_rainy_","Prob_gsg6_sd_rainy_","Prob_gsg6_q25_rainy_","Prob_gsg6_q75_rainy_","Prob_gsg6_IQR_rainy_","Prob_gsg6_ll_rainy_","Prob_gsg6_ul_rainy_"))
mainnames=as.list(c("PR gsg6_mean","PR gsg6_sd","PR gsg6_25th pct","PR gsg6_75th pct","gsg6_IQR","PR gsg6_ll","PR gsg6_ul"))
path_output <- paste0(getwd(),"/03 Output/gsg6/seasons")  
destfile<-paste0(path_output,'/')#make sure you have created a file OUTPUT in the bucket

out<-list()
for (j in 1:length(predinput))
{
  pred_val[!w] <- round(predinput[[j]],3)
  out[[j]] = setValues(mask, pred_val)
  writeRaster(out[[j]], paste0(destfile,prednames[[j]],',.tif'), overwrite=TRUE)
}
raster.list <- list.files(path=destfile,pattern =".tif$", full.names=TRUE)
raster.list <- raster.list[grep('rainy',raster.list)]
#extract raster data
rasterls<-list()
for (i in 1:length(raster.list)){
  rasterls[[i]]<-raster(raster.list[[i]])
}
d <- brick(rasterls)
plot(d)
#add informative text before the graphs
t1 = "Predicted posterior mean probability of occurrence"
t2 = "Predicted posterior lower limit probability"
t3 = "Predicted posterior upper limit probability"
t4 = "Predicted posterior inter-quartile range"
t5 = "Predicted posterior Q25"
t6 = "Predicted posterior Q75"
t7 = "Predicted posterior standard deviation"
plot.new()
pdf.options()
pdf(file=paste0(destfile,"Prediction_probability_gsg6_states_rainy.pdf"))
plot(d[[3]], main=t1)
plot(d[[2]], main=t2)
plot(d[[7]], main=t3)
plot(d[[1]], main=t4)
plot(d[[4]], main=t5)
plot(d[[5]], main=t6)
plot(d[[6]], main=t7)
dev.off();dev.off()

rb_gsg6_rainy <- d



#save outputs of the analysis (beta coefficients)

fixed_table = result_gsg6cov_rainy$summary.fixed#check that 
fixed_table <- fixed_table[,1:5]
fixed_table$sig <-ifelse(sign(fixed_table[,3])==sign(fixed_table[,5]),1,0)
fixed_table$or <- exp(fixed_table[,1])
fixed_table$or_lci <- fixed_table[,7]+(qnorm(0.025)*fixed_table[,2])
fixed_table$or_uci <- fixed_table[,7]+(qnorm(0.975)*fixed_table[,2])

write.csv(fixed_table, file="03 Output/gsg6/seasons/fixedeffects_gsg6_rainy.csv",row.names=TRUE)


########################
#### Binomial model ####
########################

#### Mesh construction
## x and y coordinates in the Response data
sample_locations = cbind(gsg6cov_cool$x, gsg6cov_cool$y)


A <- inla.spde.make.A(mesh, sample_locations) 
indexs <- inla.spde.make.index("field", spde$n.spde)



####################################################################################


covar_z = gsg6cov_cool[,-1:-8]

stk.z = inla.stack(data=list(occurence=as.vector(gsg6cov_cool$pos_gsg6), trials=as.vector(gsg6cov_cool$Nscreen)),
                   A=list(A, 1),
                   effects=list(indexs,
                                list(
                                  # z.intercept=1,
                                  z.intercept=rep(1,nrow(gsg6cov_cool)),
                                  covariate=covar_z)),
                   tag="est.z")

formula.z = as.formula(paste("occurence ~ -1 + z.intercept + f(field, model=spde) + ",paste(names(covar_z),collapse="+")))



result_gsg6cov_cool <-   inla(formula.z, family="binomial", Ntrials=inla.stack.data(stk.z)$trials,
                              data=inla.stack.data(stk.z), verbose=FALSE,
                              control.predictor = list(compute=TRUE, A=inla.stack.A(stk.z)),
                              control.fixed=list(expand.factor.strategy='inla'),
                              control.compute = list(dic=TRUE, waic=TRUE, config=TRUE),
                              control.inla = list(int.strategy="eb"))

#############################


##

pred_covar_z <- data.frame((covs_cool))
Ap <- inla.spde.make.A(mesh = mesh, loc=reference.coordinates);dim(Ap)

mask<-reference.image; NAvalue(mask)=-9999
pred_val<-getValues(mask)
w<-is.na(pred_val)
index<-1:length(w)
index<-index[!w]
pred_locs<-xyFromCell(mask,1:ncell(mask))
pred_locs<-pred_locs[!w,]
colnames(pred_locs)<-c('longitude','latitude')
locs_pred <- pred_locs


#define number of samples 
nn = 150#(1500 for final results, 150 for fast results)

#sampling the posterior
set.seed(999)
samp = inla.posterior.sample(nn, result_gsg6cov_cool)
pred = matrix(NA,nrow=dim(Ap)[1], ncol=nn)
k = dim(pred_covar_z)[2] ## number of final covariates
for (i in 1:nn){
  field = samp[[i]]$latent[grep('field',rownames(samp[[i]]$latent)),]
  intercept = samp[[i]]$latent[grep('z.intercept',rownames(samp[[i]]$latent)),]
  beta = NULL
  for (j in 1:k){
    beta[j] = samp[[i]]$latent[grep(names(pred_covar_z)[j],rownames(samp[[i]]$latent)),]
  }
  #compute beta*covariate for each covariate
  linpred<-list()
  for (j in 1:k){
    linpred[[j]]<-beta[j]*pred_covar_z[,j]
  }
  linpred<-Reduce("+",linpred)
  lp = intercept + linpred + drop(Ap%*%field)
  ## Predicted values
  pred[,i] = plogis(lp) #for a bernoulli likelihood
  # pred[,i] = exp(lp)#for a poisson likelihood
}
pred_gsg6_mean = rowMeans(pred)
pred_gsg6_ll = apply(pred, 1, function(x) quantile(x, probs=c(0.025), na.rm=TRUE))
pred_gsg6_ul = apply(pred, 1, function(x) quantile(x, probs=c(0.975), na.rm=TRUE))
pred_gsg6_sd = apply(pred, 1, sd)
pred_gsg6_25pct = apply(pred, 1, function(x) quantile(x, probs=c(0.25), na.rm=TRUE))
pred_gsg6_75pct = apply(pred, 1, function(x) quantile(x, probs=c(0.75), na.rm=TRUE))
pred_gsg6_IQR = pred_gsg6_75pct - pred_gsg6_25pct
#Saving predictive maps as raster files
predinput=list(pred_gsg6_mean,pred_gsg6_sd,pred_gsg6_25pct,pred_gsg6_75pct,pred_gsg6_IQR,pred_gsg6_ll,pred_gsg6_ul)
prednames=as.list(c("Prob_gsg6_mean_cool_","Prob_gsg6_sd_cool_","Prob_gsg6_q25_cool_","Prob_gsg6_q75_cool_","Prob_gsg6_IQR_cool_","Prob_gsg6_ll_cool_","Prob_gsg6_ul_cool_"))
mainnames=as.list(c("PR gsg6_mean","PR gsg6_sd","PR gsg6_25th pct","PR gsg6_75th pct","gsg6_IQR","PR gsg6_ll","PR gsg6_ul"))
path_output <- paste0(getwd(),"/03 Output/gsg6/seasons")  
destfile<-paste0(path_output,'/')#make sure you have created a file OUTPUT in the bucket

out<-list()
for (j in 1:length(predinput))
{
  pred_val[!w] <- round(predinput[[j]],3)
  out[[j]] = setValues(mask, pred_val)
  writeRaster(out[[j]], paste0(destfile,prednames[[j]],',.tif'), overwrite=TRUE)
}
raster.list <- list.files(path=destfile,pattern =".tif$", full.names=TRUE)
raster.list <- raster.list[grep('cool',raster.list)]
#extract raster data
rasterls<-list()
for (i in 1:length(raster.list)){
  rasterls[[i]]<-raster(raster.list[[i]])
}
d <- brick(rasterls)
plot(d)
#add informative text before the graphs
t1 = "Predicted posterior mean probability of occurrence"
t2 = "Predicted posterior lower limit probability"
t3 = "Predicted posterior upper limit probability"
t4 = "Predicted posterior inter-quartile range"
t5 = "Predicted posterior Q25"
t6 = "Predicted posterior Q75"
t7 = "Predicted posterior standard deviation"
plot.new()
pdf.options()
pdf(file=paste0(destfile,"Prediction_probability_gsg6_states_cool.pdf"))
plot(d[[3]], main=t1)
plot(d[[2]], main=t2)
plot(d[[7]], main=t3)
plot(d[[1]], main=t4)
plot(d[[4]], main=t5)
plot(d[[5]], main=t6)
plot(d[[6]], main=t7)
dev.off();dev.off()

rb_gsg6_cool <- d



#save outputs of the analysis (beta coefficients)

fixed_table = result_gsg6cov_cool$summary.fixed#check that 
fixed_table <- fixed_table[,1:5]
fixed_table$sig <-ifelse(sign(fixed_table[,3])==sign(fixed_table[,5]),1,0)
fixed_table$or <- exp(fixed_table[,1])
fixed_table$or_lci <- fixed_table[,7]+(qnorm(0.025)*fixed_table[,2])
fixed_table$or_uci <- fixed_table[,7]+(qnorm(0.975)*fixed_table[,2])

write.csv(fixed_table, file="03 Output/gsg6/seasons/fixedeffects_gsg6_cool.csv",row.names=TRUE)


