#jont modelling of sg6, csp and pcr data
#1.step install the packages you need
rm(list=ls())
#INLA used to fit Bayesian models
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

###Data

pr_dat <- read_csv("02 Cleaned data/test_sg6_csp_pcr_villages.csv")
load("02 Cleaned data/covariates_static.Rdata")
load("02 Cleaned data//pred_input_1km.Rdata")

gsg6<-data.frame(pos_gsg6=pr_dat$sg6_pos,pos_csp=pr_dat$csp_pos,pos_pcr=pr_dat$pcr_pos,Nscreen = pr_dat$n, x=pr_dat$longitude, y=pr_dat$latitude)
gsg6<-gsg6[complete.cases(gsg6),]
head(gsg6)

#**********cobine the coordinate**********************
xy <- cbind(gsg6$x, gsg6$y)

ilogit <- function(x) {1/(1+exp(-x))}

#extract covariate at point coordinates
#covariate_all<-data.frame(raster::extract(covariates_static_stack[[c(1,5,11,15,16,18,20)]], xy))
covariate_all<-data.frame(raster::extract(covariates_static_stack[[-2:-3]], xy))
covs_included <- c("accessibility_to_cities_2015_v1.0"    ,         
                   "distance_to_water",
                   "HydroSHEDS_TWI"    ,
                   "SRTM_SlopePCT_Corrected.Synoptic.Overall.Data.1km.Data"  ,
                   "tree_fraction_any_percent_2000"   ,
                   "VIIRS.DNB.v2_Clean.Background.2015.Annual.Data.1km.mean"  )
#covar_z <- data.frame(covariate_all[,covs_included])
#run VIF function to remove multicolinear variables VIF<=5
source("VIF_func.R")
keep = vif_func(data.frame(covariate_all))
## Selected covariates
covar_z = covariate_all[, keep]
#covar_z = covariate_all[, covs_included]
head(covar_z)

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


#mesh <- inla.mesh.2d(sample_locations,max.edge=c(0.1,1),cut=0.01)
#layout(1)
#plot(mesh)

spde = inla.spde2.pcmatern(mesh,alpha=2 ,
                           prior.range = c(0.03,0.1), #mean distance between villages
                           prior.sigma = c(3,0.01)
)


indexs <- inla.spde.make.index("field", spde$n.spde)
lengths(indexs)

# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
A = inla.spde.make.A(mesh=mesh, loc=as.matrix(coords));dim(A)   #A matrix

N <- 104 #number of villages



################

stack_gSG6 <- inla.stack(
  data = list(y = cbind(as.vector(gsg6$pos_gsg6), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1 <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2 <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack <- inla.stack(stack_gSG6, stack_CSP, stack_u1, stack_PCR, stack_u2) 

form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



res.z <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
              Ntrial = gsg6$Nscreen, 
              data = inla.stack.data(stack), 
              control.predictor = list(A = inla.stack.A(stack),compute=TRUE),
              control.inla = list(int.strategy = 'eb'),
              control.compute = list(config = TRUE,dic=TRUE, waic=TRUE),
              control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(res.z)



####~~~~~~~####~~~~~~~####~~~~~~~####~~~~~~~~##
#### STAGE 2: Variable selection for result ####
####~~~~~~~####~~~~~~~####~~~~~~~####~~~~~~~~##
## 1. Assess the importance of the variable one at a time

nm = noquote(names(covar_z))
DIC = rep(NA, length(nm)+1)
DIC[1] = res.z$dic$dic

nmFinal <- nm

for (i in 1:length(nm)) {
  j = nm[i]
  
  covariate_z = nmFinal[-which(nmFinal == j)]
  
  
  formula.z = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(covariate_z,collapse="+")))
  
  
  result = inla(formula.z, c("binomial","binomial","gaussian","binomial","gaussian"),
                Ntrial = gsg6$Nscreen, 
                data = inla.stack.data(stack), 
                control.predictor = list(A = inla.stack.A(stack),compute=TRUE),
                control.inla = list(int.strategy = 'eb'),
                control.compute = list(dic=TRUE, waic=TRUE, config = TRUE),
                control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))
  
  if (result$dic$dic < min(DIC[1:i])){
    nmFinal = covariate_z}
  
  DIC[i+1] = result$dic$dic
  
}


#### 2. Add variables back to the model one at a time

outnm = nm[which(!nm %in% nmFinal)]
nmFinal = nmFinal
DIC2 = rep(NA, length(outnm)+1)
DIC2[1] = min(DIC)


for (i in 1:length(outnm)) {
  j = outnm[i]
  
  covariate_z = c(nmFinal, j)
  
  formula.z = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(covariate_z,collapse="+")))
  
  result = inla(formula.z, c("binomial","binomial","gaussian","binomial","gaussian"),
                Ntrial = gsg6$Nscreen, 
                data = inla.stack.data(stack), 
                control.predictor = list(A = inla.stack.A(stack),compute=TRUE),
                control.inla = list(int.strategy = 'eb'),
                control.compute = list(dic=TRUE, waic=TRUE, config = TRUE),
                control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))
  
  if (result$dic$dic < min(DIC2[1:i])){
    nmFinal = covariate_z}
  
  DIC2[i+1] = result$dic$dic
  
}


## Compute fixed effects fromt the best model
## nmFinal is the final set of cavariates with the smallest DIC
covariate_z = nmFinal
covar_z <- covar_z[,covariate_z]
gsg6cov <- data.frame(gsg6, covar_z)
####################


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))




################

stack_gSG6 <- inla.stack(
  data = list(y = cbind(as.vector(gsg6$pos_gsg6), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1 <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2 <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack <- inla.stack(stack_gSG6, stack_CSP, stack_u1, stack_PCR, stack_u2) 

result <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
               Ntrial = gsg6$Nscreen, 
               data = inla.stack.data(stack), 
               control.predictor = list(A = inla.stack.A(stack),compute=TRUE),
               control.inla = list(int.strategy = 'eb'),
               control.compute = list(config = TRUE),
               control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(result)
## Plot results

#################

pdf("03 Output/gsg6_csp_pcr/Spatial_3outcomes.pdf")

pred_gSG6_mean <- result$summary.fitted.values[1:N,1]
pred_CSP_mean <- result$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean <- result$summary.fitted.values[((3*N)+1):(4*N),1]

png(file="03 Output/gsg6_csp_pcr/empirical_gsg6Vfitted_gsg6,empirical_cspVfitted_csp.png",    width=600, height=900)
layout(1:3)
plot((pred_gSG6_mean*100),gsg6$pos_gsg6/gsg6$Nscreen*100,xlim=c(0,100),ylim=c(0,100),xlab="Fitted Mean Prev gSG6",ylab='Observed Mean Prev gSG6')
r <- cor((pred_gSG6_mean),gsg6$pos_gsg6/gsg6$Nscreen)  
r2 <- cor((pred_gSG6_mean),gsg6$pos_gsg6/gsg6$Nscreen)^2  
mylabel = bquote(italic(r) == .(format(r, digits = 3)))
text(x = 20, y = 80, labels = mylabel)
plot((pred_CSP_mean*100),gsg6$pos_csp/gsg6$Nscreen*100,xlim=c(0,50),ylim=c(0,60),xlab="Fitted Mean Prev CSP",ylab='Observed Mean Prev CSP')
r <- cor((pred_CSP_mean),gsg6$pos_csp/gsg6$Nscreen)  
r2 <- cor((pred_CSP_mean),gsg6$pos_csp/gsg6$Nscreen)^2  
mylabel = bquote(italic(r) == .(format(r, digits = 3)))
text(x = 10, y = 50, labels = mylabel)
plot((pred_PCR_mean*100),gsg6$pos_pcr/gsg6$Nscreen*100,xlim=c(0,10),ylim=c(0,10),xlab="Fitted Mean Prev PCR",ylab='Observed Mean Prev PCR')
r <- cor((pred_PCR_mean),gsg6$pos_pcr/gsg6$Nscreen)  
r2 <- cor((pred_PCR_mean),gsg6$pos_pcr/gsg6$Nscreen)^2  
mylabel = bquote(italic(r) == .(format(r, digits = 3)))
text(x = 9, y = 8, labels = mylabel)
dev.off()


png(file="03 Output/gsg6_csp_pcr/empirical_gsg6Vcsp,fitted_gsg6Vcsp.png",    width=400, height=800)
layout(1:2)
plot(gsg6$pos_gsg6/gsg6$Nscreen*100,gsg6$pos_csp/gsg6$Nscreen*100,xlim=c(0,100),ylim=c(0,100),xlab="Observed Mean Prev gSG6",ylab='Oberved Mean Prev CSP')
r <- cor(gsg6$pos_gsg6/gsg6$Nscreen,gsg6$pos_csp/gsg6$Nscreen)  
r2 <- cor(gsg6$pos_gsg6/gsg6$Nscreen,gsg6$pos_csp/gsg6$Nscreen)^2  
mylabel = bquote(italic(r) == .(format(r, digits = 3)))
text(x = 20, y = 80, labels = mylabel)
plot((pred_gSG6_mean*100),(pred_CSP_mean*100),xlim=c(0,100),ylim=c(0,60),xlab="Fitted Mean Prev gSG6",ylab='Fitted Mean Prev CSP')
r <- cor((pred_gSG6_mean),(pred_CSP_mean))  
r2 <- cor((pred_gSG6_mean),(pred_CSP_mean))^2  
mylabel = bquote(italic(r) == .(format(r, digits = 3)))
text(x = 20, y = 40, labels = mylabel)
dev.off()

png(file="03 Output/gsg6_csp_pcr/empirical_gsg6Vpcr,fitted_gsg6Vpcr.png",    width=400, height=800)
layout(1:2)

plot(gsg6$pos_gsg6/gsg6$Nscreen*100,gsg6$pos_pcr/gsg6$Nscreen*100,xlim=c(0,100),ylim=c(0,10),xlab="Observed Mean Prev gSG6",ylab='Oberved Mean Prev PCR')
r <- cor(gsg6$pos_gsg6/gsg6$Nscreen,gsg6$pos_pcr/gsg6$Nscreen)  
r2 <- cor(gsg6$pos_gsg6/gsg6$Nscreen,gsg6$pos_pcr/gsg6$Nscreen)^2 
mylabel = bquote(italic(r) == .(format(r, digits = 3)))
text(x = 20, y = 8, labels = mylabel)
plot((pred_gSG6_mean*100),(pred_PCR_mean*100),xlim=c(0,100),ylim=c(0,10),xlab="Fitted Mean Prev gSG6",ylab='Fitted Mean Prev PCR')
r <- cor((pred_gSG6_mean),(pred_PCR_mean))  
r2 <- cor((pred_gSG6_mean),(pred_PCR_mean))^2  
mylabel = bquote(italic(r) == .(format(r, digits = 3)))
text(x = 20, y = 8, labels = mylabel)
dev.off()


png(file="03 Output/gsg6_csp_pcr/empirical_cspVpcr,fitted_cspVpcr.png",    width=400, height=800)
layout(1:2)

plot(gsg6$pos_csp/gsg6$Nscreen*100,gsg6$pos_pcr/gsg6$Nscreen*100,xlim=c(0,100),ylim=c(0,10),xlab="Observed Mean Prev CSP",ylab='Oberved Mean Prev PCR')
r <- cor(gsg6$pos_csp/gsg6$Nscreen,gsg6$pos_pcr/gsg6$Nscreen)  
r2 <- cor(gsg6$pos_csp/gsg6$Nscreen,gsg6$pos_pcr/gsg6$Nscreen)^2
mylabel = bquote(italic(r) == .(format(r, digits = 3)))
text(x = 20, y = 8, labels = mylabel)
plot((pred_CSP_mean*100),(pred_PCR_mean*100),xlim=c(0,100),ylim=c(0,10),xlab="Fitted Mean Prev CSP",ylab='Fitted Mean Prev PCR')
r <- cor((pred_CSP_mean),(pred_PCR_mean))  
r2 <- cor((pred_CSP_mean),(pred_PCR_mean))^2  
mylabel = bquote(italic(r) == .(format(r, digits = 3)))
text(x = 20, y = 8, labels = mylabel)
dev.off()




#############
#6. step: make predictions
pred_covar_z <- as.data.frame(t(covariates_static))
pred_covar_z <- pred_covar_z[,covariate_z]
#pred_covar_z <- pred_covar_z[,covs_included]
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
samp = inla.posterior.sample(nn, result)
pred = matrix(NA,nrow=dim(Ap)[1], ncol=nn)
k = dim(pred_covar_z)[2] ## number of final covariates
for (i in 1:nn){
  field = samp[[i]]$latent[grep('field_gSG6',rownames(samp[[i]]$latent)),]
  intercept = samp[[i]]$latent[grep('intercept_gSG6',rownames(samp[[i]]$latent)),]
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
prednames=as.list(c("Prob_gsg6_mean_","Prob_gsg6_sd_","Prob_gsg6_q25_","Prob_gsg6_q75_","Prob_gsg6_IQR_","Prob_gsg6_ll_","Prob_gsg6_ul_"))
mainnames=as.list(c("PR gsg6_mean","PR gsg6_sd","PR gsg6_25th pct","PR gsg6_75th pct","gsg6_IQR","PR gsg6_ll","PR gsg6_ul"))
path_output <- paste0(getwd(),"/03 Output/gsg6_csp_pcr")  
destfile<-paste0(path_output,'/')#make sure you have created a file OUTPUT in the bucket

out<-list()
for (j in 1:length(predinput))
{
  pred_val[!w] <- round(predinput[[j]],3)
  out[[j]] = setValues(mask, pred_val)
  writeRaster(out[[j]], paste0(destfile,prednames[[j]],',.tif'), overwrite=TRUE)
}
raster.list <- list.files(path=destfile,pattern =".tif$", full.names=TRUE)

raster.list <- raster.list[grep('Prob_gsg6',raster.list)]
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
pdf(file=paste0(destfile,"Prediction_probability_gsg6_states.pdf"))
plot(d[[3]], main=t1)
plot(d[[2]], main=t2)
plot(d[[7]], main=t3)
plot(d[[1]], main=t4)
plot(d[[4]], main=t5)
plot(d[[5]], main=t6)
plot(d[[6]], main=t7)
dev.off();dev.off()

rb_gsg6csppcr_gsg6 <- d


pred_gsg6 <- pred
#define number of samples 
#nn = 150#(1500 for final results, 150 for fast results)

#sampling the posterior
#set.seed(999)
#samp = inla.posterior.sample(nn, result)
pred_csp = matrix(NA,nrow=dim(Ap)[1], ncol=nn)
k = 0 ## number of final covariates
for (i in 1:nn){
  field = samp[[i]]$latent[grep('field_CSP',rownames(samp[[i]]$latent)),]
  intercept = samp[[i]]$latent[grep('intercept_CSP',rownames(samp[[i]]$latent)),]
  logit_gsg6 <- log(pred_gsg6[,i])-log(1-pred_gsg6[,i])
  #  for (j in 1:k){
  #    beta[j] = samp[[i]]$latent[grep(names(pred_covar_z)[j],rownames(samp[[i]]$latent)),]
  #  }
  #compute beta*covariate for each covariate
  #  linpred<-list()
  #  for (j in 1:k){
  #    linpred[[j]]<-beta[j]*pred_covar_z[,j]
  #  }
  #  linpred<-Reduce("+",linpred)
  lp = intercept + drop(Ap%*%field) + logit_gsg6*samp[[i]]$hyperpar["Beta for b.eta2"] #Beta for b.eta2 call in
  ## Predicted values
  pred_csp[,i] = plogis(lp) #for a bernoulli likelihood
  # pred[,i] = exp(lp)#for a poisson likelihood
}
pred_csp_mean = rowMeans(pred_csp)
pred_csp_ll = apply(pred_csp, 1, function(x) quantile(x, probs=c(0.025), na.rm=TRUE))
pred_csp_ul = apply(pred_csp, 1, function(x) quantile(x, probs=c(0.975), na.rm=TRUE))
pred_csp_sd = apply(pred_csp, 1, sd)
pred_csp_25pct = apply(pred_csp, 1, function(x) quantile(x, probs=c(0.25), na.rm=TRUE))
pred_csp_75pct = apply(pred_csp, 1, function(x) quantile(x, probs=c(0.75), na.rm=TRUE))
pred_csp_IQR = pred_csp_75pct - pred_csp_25pct
#Saving predictive maps as raster files
predinput=list(pred_csp_mean,pred_csp_sd,pred_csp_25pct,pred_csp_75pct,pred_csp_IQR,pred_csp_ll,pred_csp_ul)
prednames=as.list(c("Prob_csp_mean_","Prob_csp_sd_","Prob_csp_q25_","Prob_csp_q75_","Prob_csp_IQR_","Prob_csp_ll_","Prob_csp_ul_"))
mainnames=as.list(c("PR csp_mean","PR csp_sd","PR csp_25th pct","PR csp_75th pct","csp_IQR","PR csp_ll","PR csp_ul"))

out<-list()
for (j in 1:length(predinput))
{
  pred_val[!w] <- round(predinput[[j]],3)
  out[[j]] = setValues(mask, pred_val)
  writeRaster(out[[j]], paste0(destfile,prednames[[j]],',.tif'), overwrite=TRUE)
}
raster.list <- list.files(path=destfile,pattern =".tif$", full.names=TRUE)
raster.list <- raster.list[grep('Prob_csp',raster.list)]
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
pdf(file=paste0(destfile,"Prediction_probability_csp_states.pdf"))
plot(d[[3]], main=t1)
plot(d[[2]], main=t2)
plot(d[[7]], main=t3)
plot(d[[1]], main=t4)
plot(d[[4]], main=t5)
plot(d[[5]], main=t6)
plot(d[[6]], main=t7)
dev.off();dev.off()

rb_gsg6csppcr_csp <- d
#define number of samples 
#nn = 150#(1500 for final results, 150 for fast results)

#sampling the posterior
#set.seed(999)
#samp = inla.posterior.sample(nn, result)
pred_pcr = matrix(NA,nrow=dim(Ap)[1], ncol=nn)
k = 0 ## number of final covariates
for (i in 1:nn){
  field = samp[[i]]$latent[grep('field_PCR',rownames(samp[[i]]$latent)),]
  intercept = samp[[i]]$latent[grep('intercept_PCR',rownames(samp[[i]]$latent)),]
  logit_gsg6 <- log(pred_gsg6[,i])-log(1-pred_gsg6[,i])
  logit_csp <- log(pred_csp[,i])-log(1-pred_csp[,i])
  #  for (j in 1:k){
  #    beta[j] = samp[[i]]$latent[grep(names(pred_covar_z)[j],rownames(samp[[i]]$latent)),]
  #  }
  #compute beta*covariate for each covariate
  #  linpred<-list()
  #  for (j in 1:k){
  #    linpred[[j]]<-beta[j]*pred_covar_z[,j]
  #  }
  #  linpred<-Reduce("+",linpred)
  lp = intercept + drop(Ap%*%field) + logit_gsg6*samp[[i]]$hyperpar["Beta for b.eta3"] + #Beta for b.eta3 call in
    logit_csp*samp[[i]]$hyperpar["Beta for b.eta4"] #Beta for b.eta4 call in
  ## Predicted values
  pred_pcr[,i] = plogis(lp) #for a bernoulli likelihood
  # pred[,i] = exp(lp)#for a poisson likelihood
}
pred_pcr_mean = rowMeans(pred_pcr)
pred_pcr_ll = apply(pred_pcr, 1, function(x) quantile(x, probs=c(0.025), na.rm=TRUE))
pred_pcr_ul = apply(pred_pcr, 1, function(x) quantile(x, probs=c(0.975), na.rm=TRUE))
pred_pcr_sd = apply(pred_pcr, 1, sd)
pred_pcr_25pct = apply(pred_pcr, 1, function(x) quantile(x, probs=c(0.25), na.rm=TRUE))
pred_pcr_75pct = apply(pred_pcr, 1, function(x) quantile(x, probs=c(0.75), na.rm=TRUE))
pred_pcr_IQR = pred_pcr_75pct - pred_pcr_25pct
#Saving predictive maps as raster files
predinput=list(pred_pcr_mean,pred_pcr_sd,pred_pcr_25pct,pred_pcr_75pct,pred_pcr_IQR,pred_pcr_ll,pred_pcr_ul)
prednames=as.list(c("Prob_pcr_mean_","Prob_pcr_sd_","Prob_pcr_q25_","Prob_pcr_q75_","Prob_pcr_IQR_","Prob_pcr_ll_","Prob_pcr_ul_"))
mainnames=as.list(c("PR pcr_mean","PR pcr_sd","PR pcr_25th pct","PR pcr_75th pct","pcr_IQR","PR pcr_ll","PR pcr_ul"))

out<-list()
for (j in 1:length(predinput))
{
  pred_val[!w] <- round(predinput[[j]],3)
  out[[j]] = setValues(mask, pred_val)
  writeRaster(out[[j]], paste0(destfile,prednames[[j]],',.tif'), overwrite=TRUE)
}
raster.list <- list.files(path=destfile,pattern =".tif$", full.names=TRUE)
raster.list <- raster.list[grep('Prob_pcr',raster.list)]
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
pdf(file=paste0(destfile,"Prediction_probability_pcr_states.pdf"))
plot(d[[3]], main=t1)
plot(d[[2]], main=t2)
plot(d[[7]], main=t3)
plot(d[[1]], main=t4)
plot(d[[4]], main=t5)
plot(d[[5]], main=t6)
plot(d[[6]], main=t7)
dev.off();dev.off()

rb_gsg6csppcr_pcr <- d


save(rb_gsg6csppcr_gsg6,rb_gsg6csppcr_pcr,rb_gsg6csppcr_csp,result, file="03 Output/gsg6_csp_pcr/gsg6_csp_pcr.Rdata")
#save outputs of the analysis (beta coefficients)
tablefile<-paste0(destfile)
fixed_table = result$summary.fixed#check that 
fixed_table <- fixed_table[,1:5]
slope_table = result$summary.hyperpar[7:9,1:5]
fixed_table <- data.frame(rbind(fixed_table,slope_table))
fixed_table$sig <-ifelse(sign(fixed_table[,3])==sign(fixed_table[,5]),1,0)

fixed_table$or <- exp(fixed_table[,1])
fixed_table$or_lci <- exp(fixed_table[,1]+qnorm(0.025)*fixed_table[,2])
fixed_table$or_uci <- exp(fixed_table[,1]+qnorm(0.975)*fixed_table[,2])
write.csv(fixed_table, file=paste0("03 Output/gsg6_csp_pcr/fixedeffects_gsg6_csp_pcr_joint",".csv"),row.names=TRUE)

