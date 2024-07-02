#simple SG6 only model
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

gsg6<-data.frame(pos_gsg6=pr_dat$sg6_pos,Nscreen = pr_dat$n, x=pr_dat$longitude, y=pr_dat$latitude)
gsg6<-gsg6[complete.cases(gsg6),]
head(gsg6)

#**********cobine the coordinate**********************
xy <- cbind(gsg6$x, gsg6$y)


#extract covariate at point coordinates
covariate_all<-data.frame(raster::extract(covariates_static_stack[[-2:-3]], xy))
head(covariate_all)

#run VIF function to remove multicolinear variables VIF<=5
source("VIF_func.R")
keep = vif_func(data.frame(covariate_all))
## Selected covariates
covar_z = covariate_all[, keep]
head(covar_z)
#covar_z = covar_z[,-2] #remove access to healthcare
gsg6cov <- cbind(gsg6,covar_z)
########################
#### Binomial model ####
########################
#expand mesh to comfortably cover prediction area by adding in dummy coordinates
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


indexs <- inla.spde.make.index("field", spde$n.spde)
lengths(indexs)

# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
A = inla.spde.make.A(mesh=mesh, loc=as.matrix(coords));dim(A)   #A matrix


Y=gsg6cov$pos_gsg6

##################################
#run the model with all the covariates
#stack all data together
stk.z = inla.stack(data=list(Y=as.vector(gsg6cov$pos_gsg6), Nscreen=as.vector(gsg6cov$Nscreen)),
                   A=list(A, 1),
                   effects=list(indexs,
                                list(
                                  # z.intercept=1,
                                  z.intercept=rep(1,nrow(gsg6cov)),
                                  covariate=covar_z)),
                   tag="est.z")

#h.spec <- list(theta=list(prior='pccor1', param=c(0, 0.9)))

formula.z = as.formula(paste("Y ~ -1 + z.intercept + f(field, model=spde) + ",paste(names(covar_z),collapse="+")))

res.z = inla(formula.z, family="binomial",
             Ntrials=inla.stack.data(stk.z)$Nscreen,
             data=inla.stack.data(stk.z), verbose=F,
             control.predictor = list(compute=TRUE, A=inla.stack.A(stk.z)),
             control.fixed=list(expand.factor.strategy='inla'),
             control.compute = list(dic=TRUE, waic=TRUE, config=TRUE),
             control.inla = list(int.strategy="eb"))

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
  
  formula.z = as.formula(paste("Y ~ -1 + z.intercept + f(field, model=spde)  + ",paste(covariate_z,collapse="+")))
  
  result = inla(formula.z, family="binomial", Ntrials=inla.stack.data(stk.z)$Nscreen,
                data=inla.stack.data(stk.z), verbose=FALSE,
                control.predictor = list(compute=TRUE, A=inla.stack.A(stk.z)),
                control.fixed=list(expand.factor.strategy='inla'),
                control.compute = list(dic=TRUE, waic=TRUE, config=TRUE),
                control.inla = list(int.strategy="eb"))
  
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
  
  formula.z = as.formula(paste("Y ~ -1 + z.intercept + f(field, model=spde) + ",paste(covariate_z,collapse="+")))
  result = inla(formula.z, family="binomial", Ntrials=inla.stack.data(stk.z)$Nscreen,
                data=inla.stack.data(stk.z), verbose=FALSE,
                control.predictor = list(compute=TRUE, A=inla.stack.A(stk.z)),
                control.fixed=list(expand.factor.strategy='inla'),
                control.compute = list(dic=TRUE, waic=TRUE, config=TRUE),
                control.inla = list(int.strategy="eb"))
  
  if (result$dic$dic < min(DIC2[1:i])){
    nmFinal = covariate_z}
  
  DIC2[i+1] = result$dic$dic
  
}


## Compute fixed effects fromt the best model
## nmFinal is the final set of cavariates with the smallest DIC
covariate_z = nmFinal
covar_z <- covar_z[,covariate_z]
gsg6cov <- data.frame(gsg6, covar_z)

stk.z = inla.stack(data=list(Y=gsg6cov$pos_gsg6, Nscreen=gsg6cov$Nscreen),
                   A=list(1, A),
                   effects=list(list(
                     # z.intercept=1,
                     z.intercept=1,
                     covariate=covar_z),
                     field = indexs),
                   tag="est.z")



## Single model with Bernoulli likelihood for the occurence data
formula.z = as.formula(paste("Y ~ -1 + z.intercept +  f(field, model=spde) +",paste0(names(covar_z),collapse="+")))

result = inla(formula.z, family="binomial", Ntrials=inla.stack.data(stk.z)$Nscreen,
              data=inla.stack.data(stk.z), verbose=FALSE,
              control.predictor = list(compute=TRUE, A=inla.stack.A(stk.z)),
              control.fixed=list(expand.factor.strategy='inla'),
              control.compute = list(return.marginals=TRUE, dic=TRUE, waic=TRUE, config=TRUE),
              control.inla = list(int.strategy="eb"))
summary(result)

result$summary.hyperpar
result2<-inla.spde2.result(result, 'field', spde, do.transf=TRUE) #and then look at 
result2$summary.log.range.nominal
result2$summary.log.variance.nominal

gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
g.mean <- inla.mesh.project(gproj, result$summary.random$field$mean)
g.sd <- inla.mesh.project(gproj, result$summary.random$field$sd)
library(gridExtra)
library(lattice)
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions = heat.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd' ,col.regions = heat.colors(16)), nrow=1)



gsg6_fitted <- data.frame(cbind(Prob_gsg6_mean_.=result$summary.fitted.values$mean[1:nrow(gsg6)],
                                Prob_gsg6_sd_.=result$summary.fitted.values$sd[1:nrow(gsg6)],
                                Prob_gsg6_ll_.=result$summary.fitted.values$`0.025quant`[1:nrow(gsg6)],
                                Prob_gsg6_ul_.=result$summary.fitted.values$`0.975quant`[1:nrow(gsg6)]))

gsg6_result <- result



#############
#6. step: make predictions
pred_covar_z <- as.data.frame(t(covariates_static))
pred_covar_z <- pred_covar_z[,covariate_z]
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
prednames=as.list(c("Prob_gsg6_mean_","Prob_gsg6_sd_","Prob_gsg6_q25_","Prob_gsg6_q75_","Prob_gsg6_IQR_","Prob_gsg6_ll_","Prob_gsg6_ul_"))
mainnames=as.list(c("PR gsg6_mean","PR gsg6_sd","PR gsg6_25th pct","PR gsg6_75th pct","gsg6_IQR","PR gsg6_ll","PR gsg6_ul"))
path_output <- paste0(getwd(),"/03 Output/gsg6")  
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

gsg6_prob_mean <- 
  tm_shape(d[[3]]) +
  tm_raster( style = "fixed", breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
             palette = viridisLite::magma(10, begin = 1, end = 0.26), colorNA = NULL, title = "gSG6 %") +
  tm_layout(title = paste("gSG6 %"), 
            legend.position = c("left", "bottom"))

tmap_save(
  tm = gsg6_prob_mean,
  filename = paste0(destfile,'gsg6_prob_mean.png')
)

gsg6_prob_sd <- 
  tm_shape(d[[6]]) +
  tm_raster( style = "fixed", breaks=c(0,0.005,0.01,0.02,0.04,0.08,0.16,0.32),
             palette = viridisLite::magma(10, begin = 1, end = 0.26), colorNA = NULL, title = "gSG6 SD") +
  tm_layout(title = paste("gSG6 SD"), 
            legend.position = c("left", "bottom"))

tmap_save(
  tm = gsg6_prob_sd,
  filename = paste0(destfile,'gsg6_prob_sd.png')
)


gsg6_mean <- raster::extract(d[[3]], reference.coordinates)
gsg6_mean_prev <- gsg6_mean*100

mask<-reference.image; NAvalue(mask)=-9999
pred_val<-getValues(mask)
w<-is.na(pred_val)
pred_val[!w] <- round(gsg6_mean_prev,1)
out = setValues(mask, pred_val)

gsg6_prev_mean <- 
  tm_shape(out) +
  tm_raster( style = "fixed", breaks=c(0,10,20,30,40,50,60,70,80,90,100),
             palette = viridisLite::magma(10, begin = 1, end = 0.26), colorNA = NULL, title = "gSG6 %") +
  tm_layout(title = paste("gSG6 %"), 
            legend.position = c("left", "bottom"))

tmap_save(
  tm = gsg6_prev_mean,
  filename =paste0(destfile,'gsg6_prev_mean.png')
)
#save outputs of the analysis (beta coefficients)
tablefile<-paste0(destfile)
fixed_table = result$summary.fixed#check that 
fixed_table <- fixed_table[,1:5]
fixed_table$sig <-ifelse(sign(fixed_table[,3])==sign(fixed_table[,5]),1,0)
fixed_table$or <- exp(fixed_table[,1])
fixed_table$or_lci <- exp(fixed_table[,1]+qnorm(0.025)*fixed_table[,2])
fixed_table$or_uci <- exp(fixed_table[,1]+qnorm(0.975)*fixed_table[,2])
write.csv(fixed_table, file=paste0(tablefile,"fixedeffects_gsg6",".csv"),row.names=TRUE)

##~~##
##~~##
#1
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs1 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords1 <- xy[train_ind,]
test_coords1 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords1));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords1));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred1 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_1 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_1 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_1 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#2
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(2)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs2 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords2 <- xy[train_ind,]
test_coords2 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords2));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords2));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred2 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_2 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_2 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_2 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#3
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(3)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs3 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords3 <- xy[train_ind,]
test_coords3 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords3));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords3));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred3 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_3 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_3 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_3 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#4
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(4)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs4 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords4 <- xy[train_ind,]
test_coords4 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords4));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords4));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred4 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_4 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_4 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_4 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#5
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(5)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs5 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords5 <- xy[train_ind,]
test_coords5 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords5));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords5));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred5 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_5 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_5 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_5 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#6
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(6)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs6 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords6 <- xy[train_ind,]
test_coords6 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords6));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords6));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred6 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_6 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_6 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_6 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#7
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(7)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs7 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords7 <- xy[train_ind,]
test_coords7 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords7));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords7));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred7 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_7 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_7 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_7 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#8
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(8)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs8 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords8 <- xy[train_ind,]
test_coords8 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords8));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords8));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred8 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_8 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_8 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_8 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#9
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(9)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs9 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords9 <- xy[train_ind,]
test_coords9 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords9));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords9));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred9 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_9 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_9 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_9 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#10
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(10)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs10 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords10 <- xy[train_ind,]
test_coords10 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords10));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords10));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred10 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_10 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_10 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_10 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#########################################################
##~~##
#11
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(1123)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs11 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords11 <- xy[train_ind,]
test_coords11 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords11));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords11));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred11 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_11 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_11 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_11 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#12
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(12)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs12 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords12 <- xy[train_ind,]
test_coords12 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.12d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords12));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords12));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred12 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_12 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_12 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_12 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#13
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(13)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs13 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords13 <- xy[train_ind,]
test_coords13 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh13 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords13));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords13));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred13 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_13 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_13 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_13 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#14
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(14)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs14 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords14 <- xy[train_ind,]
test_coords14 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords14));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords14));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred14 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_14 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_14 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_14 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#15
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(15)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs15 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords15 <- xy[train_ind,]
test_coords15 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords15));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords15));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred15 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_15 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_15 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_15 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#16
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(16)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs16 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords16 <- xy[train_ind,]
test_coords16 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords16));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords16));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred16 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_16 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_16 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_16 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#17
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(17)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs17 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords17 <- xy[train_ind,]
test_coords17 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords17));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords17));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred17 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_17 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_17 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_17 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#18
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(18)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs18 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords18 <- xy[train_ind,]
test_coords18 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords18));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords18));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred18 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_18 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_18 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_18 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

#19
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(119)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs19 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords19 <- xy[train_ind,]
test_coords19 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords19));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords19));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred19 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_19 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_19 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_19 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form


#20
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.90 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(20)
train_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)

train <- gsg6cov[train_ind, ]
test <- gsg6cov[-train_ind, ]
obs20 <- cbind(test$pos_gsg6,test$Nscreen,test$pr)#this is the number pos/number examined
test$pos_gsg6 <- NA  #make the y values for test NA

#lastly, create the training and testing coordinates
train_coords20 <- xy[train_ind,]
test_coords20 <- xy[-train_ind,]


traincov<-train[,-1:-5]
testcov<-test[,-1:-5]


# mesh3 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=c(0.5,1), offset = c(0.5, 1), cutoff = 0.3)
Atrain = inla.spde.make.A(mesh=mesh, loc=as.matrix(train_coords20));dim(Atrain)   #A matrix
Atest <- inla.spde.make.A(mesh = mesh, loc=as.matrix(test_coords20));dim(Atest)


Y=gsg6cov$pos_gsg6
# stack for estimation stk.e
stk.et <- inla.stack(
  tag = "est",
  data = list(Y=train$pos_gsg6, n=train$Nscreen),
  A = list(1, Atrain),
  effects = list(list(b0 = 1, covariate=traincov), field = indexs)
)

# stack for prediction stk.p
stk.pt <- inla.stack(
  tag = "pred",
  data = list(Y = test$pos_gsg6, n = test$Nscreen),
  A = list(1, Atest),
  effects = 
    list(list(b0 = 1, covariate=testcov),
         field = indexs
    )
)

# stk.full has stk.e and stk.p
stk.fulltt <- inla.stack(stk.et, stk.pt)



#formula <- Y ~ +1 + f(spatial.field, model=spde) + Prec.z +Alt.Z + Pop.Z
formula = as.formula(paste("Y ~ -1 + b0 + f(field, model=spde) +",paste0(names(testcov),collapse="+")))
#5. step: run the model

#run model
model_gsg6_tt <- inla(formula,
                      family = "binomial", Ntrials = n,
                      control.family = list(link = "logit"),
                      control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, config = TRUE),
                      data = inla.stack.data(stk.fulltt),
                      control.predictor = list(
                        compute = TRUE, link = 1,
                        A = inla.stack.A(stk.fulltt)
                      )
)

index.pred <- inla.stack.index(stk.fulltt, "pred")$data
pred20 <- model_gsg6_tt$summary.fitted.values$mean[index.pred] #the posterior is in logit form
pred_ll_20 <- model_gsg6_tt$summary.fitted.values$`0.025quant`[index.pred] #the posterior is in logit form
pred_ul_20 <- model_gsg6_tt$summary.fitted.values$`0.975quant`[index.pred] #the posterior is in logit form
pred_sd_20 <- model_gsg6_tt$summary.fitted.values$sd[index.pred] #the posterior is in logit form

############

obs <- rbind(obs1,obs2,obs3,obs4,obs5,obs6,obs7,obs8,obs9,obs10,
             obs11,obs12,obs13,obs14,obs15,obs16,obs17,obs18,obs19,obs20)
pred <- c(pred1,pred2,pred3,pred4,pred5,pred6,pred7,pred8,pred9,pred10,
          pred11,pred12,pred13,pred14,pred15,pred16,pred17,pred18,pred19,pred20)

pred_ll <- c(pred_ll_1,pred_ll_2,pred_ll_3,pred_ll_4,pred_ll_5,pred_ll_6,pred_ll_7,pred_ll_8,pred_ll_9,pred_ll_10,
             pred_ll_11,pred_ll_12,pred_ll_13,pred_ll_14,pred_ll_15,pred_ll_16,pred_ll_17,pred_ll_18,pred_ll_19,pred_ll_20)
pred_ul <- c(pred_ul_1,pred_ul_2,pred_ul_3,pred_ul_4,pred_ul_5,pred_ul_6,pred_ul_7,pred_ul_8,pred_ul_9,pred_ul_10,
             pred_ul_11,pred_ul_12,pred_ul_13,pred_ul_14,pred_ul_15,pred_ul_16,pred_ul_17,pred_ul_18,pred_ul_19,pred_ul_20)
pred_sd <- c(pred_sd_1,pred_sd_2,pred_sd_3,pred_sd_4,pred_sd_5,pred_sd_6,pred_sd_7,pred_sd_8,pred_sd_9,pred_sd_10,
             pred_sd_11,pred_sd_12,pred_sd_13,pred_sd_14,pred_sd_15,pred_sd_16,pred_sd_17,pred_sd_18,pred_sd_19,pred_sd_20)


test_coords <- rbind(test_coords1,test_coords2,test_coords3,test_coords4,test_coords5,test_coords6,test_coords7,test_coords8,test_coords9,test_coords10,
                     test_coords11,test_coords12,test_coords13,test_coords14,test_coords15,test_coords16,test_coords17,test_coords18,test_coords19,test_coords20)

train_coords <- rbind(train_coords1,train_coords2,train_coords3,train_coords4,train_coords5,train_coords6,train_coords7,train_coords8,train_coords9,train_coords10,
                      train_coords11,train_coords12,train_coords13,train_coords14,train_coords15,train_coords16,train_coords17,train_coords18,train_coords19,train_coords20)


pred_dat <- data.frame(obs,pred,pred_sd,pred_ll,pred_ul)

obs_df <- as.data.frame(obs)

gsg6_ll_qbeta <- qbeta(0.025, shape1 = 1+obs_df$V1, shape2 = 1+(obs_df$V2-obs_df$V1))
gsg6_ul_qbeta <- qbeta(0.975, shape1 = 1+obs_df$V1, shape2 = 1+(obs_df$V2-obs_df$V1))
gsg6_med_qbeta <- qbeta(0.5, shape1 = 1+obs_df$V1, shape2 = 1+(obs_df$V2-obs_df$V1))

full_tt <- data.frame(obs=gsg6_med_qbeta, obs_ll=gsg6_ll_qbeta, obs_ul=gsg6_ul_qbeta, pred= pred, pred_ll=pred_ll, pred_ul = pred_ul)


###########################################
#################################

#we have hashed out binning by obs and now bin by pred
#pred_dat <- pred_dat %>% mutate(gsg6_bin =ntile(X3, n=10))
pred_dat <- pred_dat %>% mutate(gsg6_bin =ntile(pred, n=10))
pred_dat <- pred_dat %>% mutate(k1 =1)
pred_dat <- within(pred_dat, gsg6bin_no_villages <- ave(gsg6_bin, list(gsg6_bin, k1), FUN=length))

pred_dat <- pred_dat %>% mutate(k1 =1)
pred_dat <- within(pred_dat, gsg6bin_no_villages <- ave(gsg6_bin, list(gsg6_bin, k1), FUN=length))

head(pred_dat)

pred_dat <- pred_dat %>% 
  group_by(gsg6_bin) %>% 
  mutate(gsg6bin_sum_n = sum(X2)) %>% 
  ungroup()

pred_dat <- pred_dat %>% 
  group_by(gsg6_bin) %>% 
  mutate(gsg6bin_sum_gsg6 = sum(X1)) %>% 
  ungroup()


#################################
#gsg6 percentile weights
##calculate population weights for township
pred_dat <- pred_dat %>% mutate(pred_dat,mean_gsg6_sg6pctweight=gsg6bin_no_villages*pred*(X2/gsg6bin_sum_n))
pred_dat <- pred_dat %>% mutate(pred_dat,sd_sq_gsg6_sg6pctweight=gsg6bin_no_villages*pred_sd^2*(X2/gsg6bin_sum_n))


pred_sg6pct_weights = aggregate(pred_dat,
                                by = list(pred_dat$gsg6_bin),
                                FUN = mean)

pred_sg6pct_weights$sd_gsg6_sg6pctweight  <- sqrt(pred_sg6pct_weights$sd_sq_gsg6_sg6pctweight)


head(pred_sg6pct_weights)
##calculate 95%CredI

pred_sg6pct_weights <- pred_sg6pct_weights %>% mutate(pred_sg6pct_weights,ll_gsg6_sg6pctweight=mean_gsg6_sg6pctweight-(2*sd_gsg6_sg6pctweight))
pred_sg6pct_weights <- pred_sg6pct_weights %>% mutate(pred_sg6pct_weights,ul_gsg6_sg6pctweight=mean_gsg6_sg6pctweight+(2*sd_gsg6_sg6pctweight))


gsg6_ll_qbeta_agg <- qbeta(0.025, shape1 = 1+pred_sg6pct_weights$gsg6bin_sum_gsg6, shape2 = 1+(pred_sg6pct_weights$gsg6bin_sum_n-pred_sg6pct_weights$gsg6bin_sum_gsg6))
gsg6_ul_qbeta_agg <- qbeta(0.975, shape1 = 1+pred_sg6pct_weights$gsg6bin_sum_gsg6, shape2 = 1+(pred_sg6pct_weights$gsg6bin_sum_n-pred_sg6pct_weights$gsg6bin_sum_gsg6))
gsg6_med_qbeta_agg <- qbeta(0.5, shape1 = 1+pred_sg6pct_weights$gsg6bin_sum_gsg6, shape2 = 1+(pred_sg6pct_weights$gsg6bin_sum_n-pred_sg6pct_weights$gsg6bin_sum_gsg6))

agg_tt <- data.frame(obs=gsg6_med_qbeta_agg, obs_ll=gsg6_ll_qbeta_agg,obs_ul=gsg6_ul_qbeta_agg,
                     pred=pred_sg6pct_weights$mean_gsg6_sg6pctweight, pred_ll=pred_sg6pct_weights$ll_gsg6_sg6pctweight,pred_ul=pred_sg6pct_weights$ul_gsg6_sg6pctweight)

##adjust so no prob<0 or >1
agg_tt <- agg_tt %>% mutate(pred_ll_adj=case_when( pred_ll <0 ~ 0,
                                                   pred_ll >= 0 ~ pred_ll))


agg_tt <- agg_tt %>% mutate(pred_adj=case_when( pred >1 ~ 1,
                                                pred <0 ~ 0,
                                                pred <= 1 ~ pred))

agg_tt <- agg_tt %>% mutate(pred_ul_adj=case_when( pred_ul >1 ~ 1,
                                                   pred_ul <0 ~ 0,
                                                   pred_ul <= 1 ~ pred_ul))


ggplot() +
  geom_point(data = full_tt, mapping =aes(x = obs, y = pred), col="pink") +
  geom_pointrange(data = full_tt, mapping = aes(x = obs, y = pred, xmin = obs_ll, xmax = obs_ul), col="pink")+
  geom_pointrange(data = full_tt, mapping = aes(x = obs, y = pred, ymin = pred_ll, ymax = pred_ul), col="pink") +
  geom_point(data = agg_tt, mapping = aes(x = obs, y = pred), alpha=0.8, col="black") +
  geom_pointrange(data = agg_tt, mapping = aes(x = obs, y = pred, xmin = obs_ll, xmax = obs_ul), alpha=0.8, col="black")+
  geom_pointrange(data = agg_tt, mapping = aes(x = obs, y = pred, ymin = pred_ll_adj, ymax = pred_ul_adj), alpha=0.8, col="black") +
  scale_y_continuous("Predicted", breaks = seq(from = 0, to = 1, by = .1), limits=c(0,1)) +
  scale_x_continuous("Observed", breaks = seq(from = 0,to = 1,by = .10), limits=c(0,1)) +
  theme_classic()+ggtitle("gSG6 spatial model fit")  +
  geom_smooth(data = agg_tt, mapping =aes(x = obs, y = pred), col="black", linetype = "dashed", method = "lm", se=FALSE, fullrange=TRUE)

cor(agg_tt$obs,agg_tt$pred)

ggplot() +
  geom_point(data = full_tt, mapping =aes(x = pred, y = obs), col="pink") +
  geom_pointrange(data = full_tt, mapping = aes(x = pred, y = obs, xmin = pred_ll, xmax = pred_ul), col="pink")+
  geom_pointrange(data = full_tt, mapping = aes(x = pred, y = obs, ymin = obs_ll, ymax = obs_ul), col="pink") +
  geom_point(data = agg_tt, mapping = aes(x = pred, y = obs), alpha=0.8, col="black") +
  geom_pointrange(data = agg_tt, mapping = aes(x = pred, y = obs, xmin = pred_ll_adj, xmax = pred_ul_adj), alpha=0.8, col="black")+
  geom_pointrange(data = agg_tt, mapping = aes(x = pred, y = obs, ymin = obs_ll, ymax = obs_ul), alpha=0.8, col="black") +
  scale_y_continuous("Observed", breaks = seq(from = 0, to = 1, by = .1), limits=c(0,1)) +
  scale_x_continuous("Predicted", breaks = seq(from = 0,to = 1,by = .10), limits=c(0,1)) +
  theme_classic()+ggtitle("gSG6 spatial model fit")  +
  geom_smooth(data = agg_tt, mapping =aes(x = pred, y = obs), col="black", linetype = "dashed", method = "lm", se=FALSE, fullrange=TRUE)

cor(agg_tt$pred,agg_tt$obs)


##################################

full_tt <- full_tt %>% mutate(pred_ll_adj=case_when( pred_ll <0 ~ 0,
                                                     pred_ll >= 0 ~ pred_ll))

full_tt <- full_tt %>% mutate(pred_adj=case_when( pred >1 ~ 1,
                                                  pred <0 ~ 0,
                                                  pred <= 1 ~ pred))

full_tt <- full_tt %>% mutate(pred_ul_adj=case_when( pred_ul >1 ~ 1,
                                                     pred_ul <0 ~ 0,
                                                     pred_ul <= 1 ~ pred_ul))



agg_tt <- agg_tt*100
full_tt <- full_tt*100


ggplot() +
  geom_point(data = full_tt, mapping =aes(x = pred_adj, y = obs), col="pink") +
  geom_pointrange(data = full_tt, mapping = aes(x = pred_adj, y = obs, xmin = pred_ll, xmax = pred_ul), col="pink")+
  geom_pointrange(data = full_tt, mapping = aes(x = pred_adj, y = obs, ymin = obs_ll, ymax = obs_ul), col="pink") +
  geom_point(data = agg_tt, mapping = aes(x = pred_adj, y = obs), alpha=0.8, col="black") +
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_adj, y = obs, xmin = pred_ll_adj, xmax = pred_ul_adj), alpha=0.8, col="black")+
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_adj, y = obs, ymin = obs_ll, ymax = obs_ul), alpha=0.8, col="black") +
  scale_y_continuous("Observed gSG6 Prevalence", breaks = seq(from = 0, to = 100, by = 20), limits=c(0,100)) +
  scale_x_continuous("Predicted gSG6 Prevalence", breaks = seq(from = 0,to = 100,by = 20), limits=c(0,100)) +
  theme_classic()+ggtitle("gSG6 spatial model fit")  +
  geom_smooth(data = agg_tt, mapping =aes(x = pred_adj, y = obs), col="black", linetype = "dashed", method = "lm", se=FALSE, fullrange=TRUE)

cor(agg_tt$obs,agg_tt$pred_adj)
cor(full_tt$obs,full_tt$pred)

