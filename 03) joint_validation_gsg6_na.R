
library(naniar)

#1
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(1)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df1 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df1$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt1 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                       Ntrial = gsg6cov$Nscreen, 
                       data = inla.stack.data(stack_test), 
                       control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                       control.inla = list(int.strategy = 'eb'),
                       control.compute = list(config = TRUE),
                       control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt1)


pred_gSG6_mean_tt1 <- joint_gsg6_tt1$summary.fitted.values[1:N,1]
pred_CSP_mean_tt1 <- joint_gsg6_tt1$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt1 <- joint_gsg6_tt1$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt1 <- joint_gsg6_tt1$summary.fitted.values[1:N,2]
pred_CSP_sd_tt1 <- joint_gsg6_tt1$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt1 <- joint_gsg6_tt1$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt1 <- joint_gsg6_tt1$summary.fitted.values[1:N,3]
pred_CSP_ll_tt1 <- joint_gsg6_tt1$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt1 <- joint_gsg6_tt1$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt1 <- joint_gsg6_tt1$summary.fitted.values[1:N,5]
pred_CSP_ul_tt1 <- joint_gsg6_tt1$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt1 <- joint_gsg6_tt1$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt1 <- cbind(pred_gSG6_mean_tt1, pred_gSG6_sd_tt1, pred_gSG6_ll_tt1, pred_gSG6_ul_tt1, 
                  pred_CSP_mean_tt1 , pred_CSP_sd_tt1, pred_CSP_ll_tt1 , pred_CSP_ul_tt1, 
                  pred_PCR_mean_tt1, pred_PCR_sd_tt1,  pred_PCR_ll_tt1, pred_PCR_ul_tt1)

pred_tt1 <-pred_tt1[test_ind,]
obs_tt1 <- gsg6cov[test_ind,c(1:4,7)]


#2
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(2)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df2 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df2$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt2 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                       Ntrial = gsg6cov$Nscreen, 
                       data = inla.stack.data(stack_test), 
                       control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                       control.inla = list(int.strategy = 'eb'),
                       control.compute = list(config = TRUE),
                       control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt2)


pred_gSG6_mean_tt2 <- joint_gsg6_tt2$summary.fitted.values[1:N,1]
pred_CSP_mean_tt2 <- joint_gsg6_tt2$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt2 <- joint_gsg6_tt2$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt2 <- joint_gsg6_tt2$summary.fitted.values[1:N,2]
pred_CSP_sd_tt2 <- joint_gsg6_tt2$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt2 <- joint_gsg6_tt2$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt2 <- joint_gsg6_tt2$summary.fitted.values[1:N,3]
pred_CSP_ll_tt2 <- joint_gsg6_tt2$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt2 <- joint_gsg6_tt2$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt2 <- joint_gsg6_tt2$summary.fitted.values[1:N,5]
pred_CSP_ul_tt2 <- joint_gsg6_tt2$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt2 <- joint_gsg6_tt2$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt2 <- cbind(pred_gSG6_mean_tt2, pred_gSG6_sd_tt2, pred_gSG6_ll_tt2, pred_gSG6_ul_tt2, 
                  pred_CSP_mean_tt2 , pred_CSP_sd_tt2, pred_CSP_ll_tt2 , pred_CSP_ul_tt2, 
                  pred_PCR_mean_tt2, pred_PCR_sd_tt2,  pred_PCR_ll_tt2, pred_PCR_ul_tt2)

pred_tt2 <-pred_tt2[test_ind,]
obs_tt2 <- gsg6cov[test_ind,c(1:4,7)]

#3
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(3)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df3 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df3$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt3 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                       Ntrial = gsg6cov$Nscreen, 
                       data = inla.stack.data(stack_test), 
                       control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                       control.inla = list(int.strategy = 'eb'),
                       control.compute = list(config = TRUE),
                       control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt3)


pred_gSG6_mean_tt3 <- joint_gsg6_tt3$summary.fitted.values[1:N,1]
pred_CSP_mean_tt3 <- joint_gsg6_tt3$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt3 <- joint_gsg6_tt3$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt3 <- joint_gsg6_tt3$summary.fitted.values[1:N,2]
pred_CSP_sd_tt3 <- joint_gsg6_tt3$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt3 <- joint_gsg6_tt3$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt3 <- joint_gsg6_tt3$summary.fitted.values[1:N,3]
pred_CSP_ll_tt3 <- joint_gsg6_tt3$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt3 <- joint_gsg6_tt3$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt3 <- joint_gsg6_tt3$summary.fitted.values[1:N,5]
pred_CSP_ul_tt3 <- joint_gsg6_tt3$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt3 <- joint_gsg6_tt3$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt3 <- cbind(pred_gSG6_mean_tt3, pred_gSG6_sd_tt3, pred_gSG6_ll_tt3, pred_gSG6_ul_tt3, 
                  pred_CSP_mean_tt3 , pred_CSP_sd_tt3, pred_CSP_ll_tt3 , pred_CSP_ul_tt3, 
                  pred_PCR_mean_tt3, pred_PCR_sd_tt3,  pred_PCR_ll_tt3, pred_PCR_ul_tt3)

pred_tt3 <-pred_tt3[test_ind,]
obs_tt3 <- gsg6cov[test_ind,c(1:4,7)]

#4
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(4)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df4 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df4$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt4 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                       Ntrial = gsg6cov$Nscreen, 
                       data = inla.stack.data(stack_test), 
                       control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                       control.inla = list(int.strategy = 'eb'),
                       control.compute = list(config = TRUE),
                       control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt4)


pred_gSG6_mean_tt4 <- joint_gsg6_tt4$summary.fitted.values[1:N,1]
pred_CSP_mean_tt4 <- joint_gsg6_tt4$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt4 <- joint_gsg6_tt4$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt4 <- joint_gsg6_tt4$summary.fitted.values[1:N,2]
pred_CSP_sd_tt4 <- joint_gsg6_tt4$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt4 <- joint_gsg6_tt4$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt4 <- joint_gsg6_tt4$summary.fitted.values[1:N,3]
pred_CSP_ll_tt4 <- joint_gsg6_tt4$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt4 <- joint_gsg6_tt4$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt4 <- joint_gsg6_tt4$summary.fitted.values[1:N,5]
pred_CSP_ul_tt4 <- joint_gsg6_tt4$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt4 <- joint_gsg6_tt4$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt4 <- cbind(pred_gSG6_mean_tt4, pred_gSG6_sd_tt4, pred_gSG6_ll_tt4, pred_gSG6_ul_tt4, 
                  pred_CSP_mean_tt4 , pred_CSP_sd_tt4, pred_CSP_ll_tt4 , pred_CSP_ul_tt4, 
                  pred_PCR_mean_tt4, pred_PCR_sd_tt4,  pred_PCR_ll_tt4, pred_PCR_ul_tt4)

pred_tt4 <-pred_tt4[test_ind,]
obs_tt4 <- gsg6cov[test_ind,c(1:4,7)]
#5
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(5)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df5 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df5$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt5 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                       Ntrial = gsg6cov$Nscreen, 
                       data = inla.stack.data(stack_test), 
                       control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                       control.inla = list(int.strategy = 'eb'),
                       control.compute = list(config = TRUE),
                       control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt5)


pred_gSG6_mean_tt5 <- joint_gsg6_tt5$summary.fitted.values[1:N,1]
pred_CSP_mean_tt5 <- joint_gsg6_tt5$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt5 <- joint_gsg6_tt5$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt5 <- joint_gsg6_tt5$summary.fitted.values[1:N,2]
pred_CSP_sd_tt5 <- joint_gsg6_tt5$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt5 <- joint_gsg6_tt5$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt5 <- joint_gsg6_tt5$summary.fitted.values[1:N,3]
pred_CSP_ll_tt5 <- joint_gsg6_tt5$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt5 <- joint_gsg6_tt5$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt5 <- joint_gsg6_tt5$summary.fitted.values[1:N,5]
pred_CSP_ul_tt5 <- joint_gsg6_tt5$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt5 <- joint_gsg6_tt5$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt5 <- cbind(pred_gSG6_mean_tt5, pred_gSG6_sd_tt5, pred_gSG6_ll_tt5, pred_gSG6_ul_tt5, 
                  pred_CSP_mean_tt5 , pred_CSP_sd_tt5, pred_CSP_ll_tt5 , pred_CSP_ul_tt5, 
                  pred_PCR_mean_tt5, pred_PCR_sd_tt5,  pred_PCR_ll_tt5, pred_PCR_ul_tt5)

pred_tt5 <-pred_tt5[test_ind,]
obs_tt5 <- gsg6cov[test_ind,c(1:4,7)]

#6
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(6)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df6 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df6$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt6 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                       Ntrial = gsg6cov$Nscreen, 
                       data = inla.stack.data(stack_test), 
                       control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                       control.inla = list(int.strategy = 'eb'),
                       control.compute = list(config = TRUE),
                       control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt6)


pred_gSG6_mean_tt6 <- joint_gsg6_tt6$summary.fitted.values[1:N,1]
pred_CSP_mean_tt6 <- joint_gsg6_tt6$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt6 <- joint_gsg6_tt6$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt6 <- joint_gsg6_tt6$summary.fitted.values[1:N,2]
pred_CSP_sd_tt6 <- joint_gsg6_tt6$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt6 <- joint_gsg6_tt6$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt6 <- joint_gsg6_tt6$summary.fitted.values[1:N,3]
pred_CSP_ll_tt6 <- joint_gsg6_tt6$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt6 <- joint_gsg6_tt6$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt6 <- joint_gsg6_tt6$summary.fitted.values[1:N,5]
pred_CSP_ul_tt6 <- joint_gsg6_tt6$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt6 <- joint_gsg6_tt6$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt6 <- cbind(pred_gSG6_mean_tt6, pred_gSG6_sd_tt6, pred_gSG6_ll_tt6, pred_gSG6_ul_tt6, 
                  pred_CSP_mean_tt6 , pred_CSP_sd_tt6, pred_CSP_ll_tt6 , pred_CSP_ul_tt6, 
                  pred_PCR_mean_tt6, pred_PCR_sd_tt6,  pred_PCR_ll_tt6, pred_PCR_ul_tt6)

pred_tt6 <-pred_tt6[test_ind,]
obs_tt6 <- gsg6cov[test_ind,c(1:4,7)]
#7
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(7)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df7 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df7$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt7 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                       Ntrial = gsg6cov$Nscreen, 
                       data = inla.stack.data(stack_test), 
                       control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                       control.inla = list(int.strategy = 'eb'),
                       control.compute = list(config = TRUE),
                       control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt7)


pred_gSG6_mean_tt7 <- joint_gsg6_tt7$summary.fitted.values[1:N,1]
pred_CSP_mean_tt7 <- joint_gsg6_tt7$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt7 <- joint_gsg6_tt7$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt7 <- joint_gsg6_tt7$summary.fitted.values[1:N,2]
pred_CSP_sd_tt7 <- joint_gsg6_tt7$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt7 <- joint_gsg6_tt7$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt7 <- joint_gsg6_tt7$summary.fitted.values[1:N,3]
pred_CSP_ll_tt7 <- joint_gsg6_tt7$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt7 <- joint_gsg6_tt7$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt7 <- joint_gsg6_tt7$summary.fitted.values[1:N,5]
pred_CSP_ul_tt7 <- joint_gsg6_tt7$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt7 <- joint_gsg6_tt7$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt7 <- cbind(pred_gSG6_mean_tt7, pred_gSG6_sd_tt7, pred_gSG6_ll_tt7, pred_gSG6_ul_tt7, 
                  pred_CSP_mean_tt7 , pred_CSP_sd_tt7, pred_CSP_ll_tt7 , pred_CSP_ul_tt7, 
                  pred_PCR_mean_tt7, pred_PCR_sd_tt7,  pred_PCR_ll_tt7, pred_PCR_ul_tt7)

pred_tt7 <-pred_tt7[test_ind,]
obs_tt7 <- gsg6cov[test_ind,c(1:4,7)]

#8
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(8)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df8 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df8$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt8 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                       Ntrial = gsg6cov$Nscreen, 
                       data = inla.stack.data(stack_test), 
                       control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                       control.inla = list(int.strategy = 'eb'),
                       control.compute = list(config = TRUE),
                       control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt8)


pred_gSG6_mean_tt8 <- joint_gsg6_tt8$summary.fitted.values[1:N,1]
pred_CSP_mean_tt8 <- joint_gsg6_tt8$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt8 <- joint_gsg6_tt8$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt8 <- joint_gsg6_tt8$summary.fitted.values[1:N,2]
pred_CSP_sd_tt8 <- joint_gsg6_tt8$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt8 <- joint_gsg6_tt8$summary.fitted.values[((3*N)+1):(4*N),2]

pred_gSG6_ll_tt8 <- joint_gsg6_tt8$summary.fitted.values[1:N,3]
pred_CSP_ll_tt8 <- joint_gsg6_tt8$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt8 <- joint_gsg6_tt8$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt8 <- joint_gsg6_tt8$summary.fitted.values[1:N,5]
pred_CSP_ul_tt8 <- joint_gsg6_tt8$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt8 <- joint_gsg6_tt8$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt8 <- cbind(pred_gSG6_mean_tt8, pred_gSG6_sd_tt8, pred_gSG6_ll_tt8, pred_gSG6_ul_tt8, 
                  pred_CSP_mean_tt8 , pred_CSP_sd_tt8, pred_CSP_ll_tt8 , pred_CSP_ul_tt8, 
                  pred_PCR_mean_tt8, pred_PCR_sd_tt8,  pred_PCR_ll_tt8, pred_PCR_ul_tt8)

pred_tt8 <-pred_tt8[test_ind,]
obs_tt8 <- gsg6cov[test_ind,c(1:4,7)]

#9
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(9)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df9 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df9$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt9 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                       Ntrial = gsg6cov$Nscreen, 
                       data = inla.stack.data(stack_test), 
                       control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                       control.inla = list(int.strategy = 'eb'),
                       control.compute = list(config = TRUE),
                       control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt9)


pred_gSG6_mean_tt9 <- joint_gsg6_tt9$summary.fitted.values[1:N,1]
pred_CSP_mean_tt9 <- joint_gsg6_tt9$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt9 <- joint_gsg6_tt9$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt9 <- joint_gsg6_tt9$summary.fitted.values[1:N,2]
pred_CSP_sd_tt9 <- joint_gsg6_tt9$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt9 <- joint_gsg6_tt9$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt9 <- joint_gsg6_tt9$summary.fitted.values[1:N,3]
pred_CSP_ll_tt9 <- joint_gsg6_tt9$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt9 <- joint_gsg6_tt9$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt9 <- joint_gsg6_tt9$summary.fitted.values[1:N,5]
pred_CSP_ul_tt9 <- joint_gsg6_tt9$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt9 <- joint_gsg6_tt9$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt9 <- cbind(pred_gSG6_mean_tt9, pred_gSG6_sd_tt9, pred_gSG6_ll_tt9, pred_gSG6_ul_tt9, 
                  pred_CSP_mean_tt9 , pred_CSP_sd_tt9, pred_CSP_ll_tt9 , pred_CSP_ul_tt9, 
                  pred_PCR_mean_tt9, pred_PCR_sd_tt9,  pred_PCR_ll_tt9, pred_PCR_ul_tt9)
pred_tt9 <-pred_tt9[test_ind,]
obs_tt9 <- gsg6cov[test_ind,c(1:4,7)]

#10
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(10)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df10 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df10$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt10 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt10)


pred_gSG6_mean_tt10 <- joint_gsg6_tt10$summary.fitted.values[1:N,1]
pred_CSP_mean_tt10 <- joint_gsg6_tt10$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt10 <- joint_gsg6_tt10$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt10 <- joint_gsg6_tt10$summary.fitted.values[1:N,2]
pred_CSP_sd_tt10 <- joint_gsg6_tt10$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt10 <- joint_gsg6_tt10$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt10 <- joint_gsg6_tt10$summary.fitted.values[1:N,3]
pred_CSP_ll_tt10 <- joint_gsg6_tt10$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt10 <- joint_gsg6_tt10$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt10 <- joint_gsg6_tt10$summary.fitted.values[1:N,5]
pred_CSP_ul_tt10 <- joint_gsg6_tt10$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt10 <- joint_gsg6_tt10$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt10 <- cbind(pred_gSG6_mean_tt10, pred_gSG6_sd_tt10, pred_gSG6_ll_tt10, pred_gSG6_ul_tt10, 
                   pred_CSP_mean_tt10 , pred_CSP_sd_tt10, pred_CSP_ll_tt10 , pred_CSP_ul_tt10, 
                   pred_PCR_mean_tt10, pred_PCR_sd_tt10,  pred_PCR_ll_tt10, pred_PCR_ul_tt10)

pred_tt10 <-pred_tt10[test_ind,]
obs_tt10 <- gsg6cov[test_ind,c(1:4,7)]
#11
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(11)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df11 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df11$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt11 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt11)


pred_gSG6_mean_tt11 <- joint_gsg6_tt11$summary.fitted.values[1:N,1]
pred_CSP_mean_tt11 <- joint_gsg6_tt11$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt11 <- joint_gsg6_tt11$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt11 <- joint_gsg6_tt11$summary.fitted.values[1:N,2]
pred_CSP_sd_tt11 <- joint_gsg6_tt11$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt11 <- joint_gsg6_tt11$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt11 <- joint_gsg6_tt11$summary.fitted.values[1:N,3]
pred_CSP_ll_tt11 <- joint_gsg6_tt11$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt11 <- joint_gsg6_tt11$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt11 <- joint_gsg6_tt11$summary.fitted.values[1:N,5]
pred_CSP_ul_tt11 <- joint_gsg6_tt11$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt11 <- joint_gsg6_tt11$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt11 <- cbind(pred_gSG6_mean_tt11, pred_gSG6_sd_tt11, pred_gSG6_ll_tt11, pred_gSG6_ul_tt11, 
                   pred_CSP_mean_tt11 , pred_CSP_sd_tt11, pred_CSP_ll_tt11 , pred_CSP_ul_tt11, 
                   pred_PCR_mean_tt11, pred_PCR_sd_tt11,  pred_PCR_ll_tt11, pred_PCR_ul_tt11)

pred_tt11 <-pred_tt11[test_ind,]
obs_tt11 <- gsg6cov[test_ind,c(1:4,7)]
#12
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(12)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df12 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df12$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt12 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt12)


pred_gSG6_mean_tt12 <- joint_gsg6_tt12$summary.fitted.values[1:N,1]
pred_CSP_mean_tt12 <- joint_gsg6_tt12$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt12 <- joint_gsg6_tt12$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt12 <- joint_gsg6_tt12$summary.fitted.values[1:N,2]
pred_CSP_sd_tt12 <- joint_gsg6_tt12$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt12 <- joint_gsg6_tt12$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt12 <- joint_gsg6_tt12$summary.fitted.values[1:N,3]
pred_CSP_ll_tt12 <- joint_gsg6_tt12$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt12 <- joint_gsg6_tt12$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt12 <- joint_gsg6_tt12$summary.fitted.values[1:N,5]
pred_CSP_ul_tt12 <- joint_gsg6_tt12$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt12 <- joint_gsg6_tt12$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt12 <- cbind(pred_gSG6_mean_tt12, pred_gSG6_sd_tt12, pred_gSG6_ll_tt12, pred_gSG6_ul_tt12, 
                   pred_CSP_mean_tt12 , pred_CSP_sd_tt12, pred_CSP_ll_tt12 , pred_CSP_ul_tt12, 
                   pred_PCR_mean_tt12, pred_PCR_sd_tt12,  pred_PCR_ll_tt12, pred_PCR_ul_tt12)

pred_tt12 <-pred_tt12[test_ind,]
obs_tt12 <- gsg6cov[test_ind,c(1:4,7)]
#13
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(13)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df13 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df13$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt13 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt13)


pred_gSG6_mean_tt13 <- joint_gsg6_tt13$summary.fitted.values[1:N,1]
pred_CSP_mean_tt13 <- joint_gsg6_tt13$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt13 <- joint_gsg6_tt13$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt13 <- joint_gsg6_tt13$summary.fitted.values[1:N,2]
pred_CSP_sd_tt13 <- joint_gsg6_tt13$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt13 <- joint_gsg6_tt13$summary.fitted.values[((3*N)+1):(4*N),2]

pred_gSG6_ll_tt13 <- joint_gsg6_tt13$summary.fitted.values[1:N,3]
pred_CSP_ll_tt13 <- joint_gsg6_tt13$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt13 <- joint_gsg6_tt13$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt13 <- joint_gsg6_tt13$summary.fitted.values[1:N,5]
pred_CSP_ul_tt13 <- joint_gsg6_tt13$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt13 <- joint_gsg6_tt13$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt13 <- cbind(pred_gSG6_mean_tt13, pred_gSG6_sd_tt13, pred_gSG6_ll_tt13, pred_gSG6_ul_tt13, 
                   pred_CSP_mean_tt13 , pred_CSP_sd_tt13, pred_CSP_ll_tt13 , pred_CSP_ul_tt13, 
                   pred_PCR_mean_tt13, pred_PCR_sd_tt13,  pred_PCR_ll_tt13, pred_PCR_ul_tt13)

pred_tt13 <-pred_tt13[test_ind,]
obs_tt13 <- gsg6cov[test_ind,c(1:4,7)]
#14
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(14)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df14 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df14$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt14 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt14)


pred_gSG6_mean_tt14 <- joint_gsg6_tt14$summary.fitted.values[1:N,1]
pred_CSP_mean_tt14 <- joint_gsg6_tt14$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt14 <- joint_gsg6_tt14$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt14 <- joint_gsg6_tt14$summary.fitted.values[1:N,2]
pred_CSP_sd_tt14 <- joint_gsg6_tt14$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt14 <- joint_gsg6_tt14$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt14 <- joint_gsg6_tt14$summary.fitted.values[1:N,3]
pred_CSP_ll_tt14 <- joint_gsg6_tt14$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt14 <- joint_gsg6_tt14$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt14 <- joint_gsg6_tt14$summary.fitted.values[1:N,5]
pred_CSP_ul_tt14 <- joint_gsg6_tt14$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt14 <- joint_gsg6_tt14$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt14 <- cbind(pred_gSG6_mean_tt14, pred_gSG6_sd_tt14, pred_gSG6_ll_tt14, pred_gSG6_ul_tt14, 
                   pred_CSP_mean_tt14 , pred_CSP_sd_tt14, pred_CSP_ll_tt14 , pred_CSP_ul_tt14, 
                   pred_PCR_mean_tt14, pred_PCR_sd_tt14,  pred_PCR_ll_tt14, pred_PCR_ul_tt14)

pred_tt14 <-pred_tt14[test_ind,]
obs_tt14 <- gsg6cov[test_ind,c(1:4,7)]
#15
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(15)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df15 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df15$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt15 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt15)


pred_gSG6_mean_tt15 <- joint_gsg6_tt15$summary.fitted.values[1:N,1]
pred_CSP_mean_tt15 <- joint_gsg6_tt15$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt15 <- joint_gsg6_tt15$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt15 <- joint_gsg6_tt15$summary.fitted.values[1:N,2]
pred_CSP_sd_tt15 <- joint_gsg6_tt15$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt15 <- joint_gsg6_tt15$summary.fitted.values[((3*N)+1):(4*N),2]

pred_gSG6_ll_tt15 <- joint_gsg6_tt15$summary.fitted.values[1:N,3]
pred_CSP_ll_tt15 <- joint_gsg6_tt15$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt15 <- joint_gsg6_tt15$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt15 <- joint_gsg6_tt15$summary.fitted.values[1:N,5]
pred_CSP_ul_tt15 <- joint_gsg6_tt15$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt15 <- joint_gsg6_tt15$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt15 <- cbind(pred_gSG6_mean_tt15, pred_gSG6_sd_tt15, pred_gSG6_ll_tt15, pred_gSG6_ul_tt15, 
                   pred_CSP_mean_tt15 , pred_CSP_sd_tt15, pred_CSP_ll_tt15 , pred_CSP_ul_tt15, 
                   pred_PCR_mean_tt15, pred_PCR_sd_tt15,  pred_PCR_ll_tt15, pred_PCR_ul_tt15)

pred_tt15 <-pred_tt15[test_ind,]
obs_tt15 <- gsg6cov[test_ind,c(1:4,7)]

#16
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(16)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df16 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df16$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt16 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt16)


pred_gSG6_mean_tt16 <- joint_gsg6_tt16$summary.fitted.values[1:N,1]
pred_CSP_mean_tt16 <- joint_gsg6_tt16$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt16 <- joint_gsg6_tt16$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt16 <- joint_gsg6_tt16$summary.fitted.values[1:N,2]
pred_CSP_sd_tt16 <- joint_gsg6_tt16$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt16 <- joint_gsg6_tt16$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt16 <- joint_gsg6_tt16$summary.fitted.values[1:N,3]
pred_CSP_ll_tt16 <- joint_gsg6_tt16$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt16 <- joint_gsg6_tt16$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt16 <- joint_gsg6_tt16$summary.fitted.values[1:N,5]
pred_CSP_ul_tt16 <- joint_gsg6_tt16$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt16 <- joint_gsg6_tt16$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt16 <- cbind(pred_gSG6_mean_tt16, pred_gSG6_sd_tt16, pred_gSG6_ll_tt16, pred_gSG6_ul_tt16, 
                   pred_CSP_mean_tt16 , pred_CSP_sd_tt16, pred_CSP_ll_tt16 , pred_CSP_ul_tt16, 
                   pred_PCR_mean_tt16, pred_PCR_sd_tt16,  pred_PCR_ll_tt16, pred_PCR_ul_tt16)
pred_tt16 <-pred_tt16[test_ind,]
obs_tt16 <- gsg6cov[test_ind,c(1:4,7)]

#17
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(17)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df17 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df17$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt17 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt17)


pred_gSG6_mean_tt17 <- joint_gsg6_tt17$summary.fitted.values[1:N,1]
pred_CSP_mean_tt17 <- joint_gsg6_tt17$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt17 <- joint_gsg6_tt17$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt17 <- joint_gsg6_tt17$summary.fitted.values[1:N,2]
pred_CSP_sd_tt17 <- joint_gsg6_tt17$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt17 <- joint_gsg6_tt17$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt17 <- joint_gsg6_tt17$summary.fitted.values[1:N,3]
pred_CSP_ll_tt17 <- joint_gsg6_tt17$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt17 <- joint_gsg6_tt17$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt17 <- joint_gsg6_tt17$summary.fitted.values[1:N,5]
pred_CSP_ul_tt17 <- joint_gsg6_tt17$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt17 <- joint_gsg6_tt17$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt17 <- cbind(pred_gSG6_mean_tt17, pred_gSG6_sd_tt17, pred_gSG6_ll_tt17, pred_gSG6_ul_tt17, 
                   pred_CSP_mean_tt17 , pred_CSP_sd_tt17, pred_CSP_ll_tt17 , pred_CSP_ul_tt17, 
                   pred_PCR_mean_tt17, pred_PCR_sd_tt17,  pred_PCR_ll_tt17, pred_PCR_ul_tt17)
pred_tt17 <-pred_tt17[test_ind,]

obs_tt17 <- gsg6cov[test_ind,c(1:4,7)]

#18
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(18)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df18 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df18$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt18 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt18)


pred_gSG6_mean_tt18 <- joint_gsg6_tt18$summary.fitted.values[1:N,1]
pred_CSP_mean_tt18 <- joint_gsg6_tt18$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt18 <- joint_gsg6_tt18$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt18 <- joint_gsg6_tt18$summary.fitted.values[1:N,2]
pred_CSP_sd_tt18 <- joint_gsg6_tt18$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt18 <- joint_gsg6_tt18$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt18 <- joint_gsg6_tt18$summary.fitted.values[1:N,3]
pred_CSP_ll_tt18 <- joint_gsg6_tt18$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt18 <- joint_gsg6_tt18$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt18 <- joint_gsg6_tt18$summary.fitted.values[1:N,5]
pred_CSP_ul_tt18 <- joint_gsg6_tt18$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt18 <- joint_gsg6_tt18$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt18 <- cbind(pred_gSG6_mean_tt18, pred_gSG6_sd_tt18, pred_gSG6_ll_tt18, pred_gSG6_ul_tt18, 
                   pred_CSP_mean_tt18 , pred_CSP_sd_tt18, pred_CSP_ll_tt18 , pred_CSP_ul_tt18, 
                   pred_PCR_mean_tt18, pred_PCR_sd_tt18,  pred_PCR_ll_tt18, pred_PCR_ul_tt18)

pred_tt18 <-pred_tt18[test_ind,]
obs_tt18 <- gsg6cov[test_ind,c(1:4,7)]

#19
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(19)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df19 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df19$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt19 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt19)


pred_gSG6_mean_tt19 <- joint_gsg6_tt19$summary.fitted.values[1:N,1]
pred_CSP_mean_tt19 <- joint_gsg6_tt19$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt19 <- joint_gsg6_tt19$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt19 <- joint_gsg6_tt19$summary.fitted.values[1:N,2]
pred_CSP_sd_tt19 <- joint_gsg6_tt19$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt19 <- joint_gsg6_tt19$summary.fitted.values[((3*N)+1):(4*N),2]


pred_gSG6_ll_tt19 <- joint_gsg6_tt19$summary.fitted.values[1:N,3]
pred_CSP_ll_tt19 <- joint_gsg6_tt19$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt19 <- joint_gsg6_tt19$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt19 <- joint_gsg6_tt19$summary.fitted.values[1:N,5]
pred_CSP_ul_tt19 <- joint_gsg6_tt19$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt19 <- joint_gsg6_tt19$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt19 <- cbind(pred_gSG6_mean_tt19, pred_gSG6_sd_tt19, pred_gSG6_ll_tt19, pred_gSG6_ul_tt19, 
                   pred_CSP_mean_tt19 , pred_CSP_sd_tt19, pred_CSP_ll_tt19 , pred_CSP_ul_tt19, 
                   pred_PCR_mean_tt19, pred_PCR_sd_tt19,  pred_PCR_ll_tt19, pred_PCR_ul_tt19)

pred_tt19 <-pred_tt19[test_ind,]
obs_tt19 <- gsg6cov[test_ind,c(1:4,7)]

#20
gsg6cov <- data.frame(gsg6, pr=gsg6$pos_gsg6/gsg6$Nscreen, covar_z) 
##test and train
## 75% of the sample size
smp_size <- floor(0.10 * nrow(gsg6cov))

## set the seed to make your partition reproducible
set.seed(20)
test_ind <- sample(seq_len(nrow(gsg6cov)), size = smp_size, replace = FALSE)



tt_df20 <-  gsg6cov %>% replace_with_na(replace = list(pos_gsg6 = test_ind) )

gsg6cov <- cbind(gsg6cov,pos_gsg6_test=tt_df20$pos_gsg6)

############


form = as.formula(paste("y ~ 0 + intercept_gSG6 + intercept_CSP + intercept_PCR + 
  f(field_gSG6, model = spde) + f(field_CSP, model = spde) + f(field_PCR, model = spde) +
  f(u1, w1, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(b.eta2, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(u2, w2, model='iid', hyper = list(prec = list(initial = -6, fixed=TRUE))) +
  f(b.eta3, copy='u1', hyper = list(beta = list(fixed = FALSE))) +
  f(b.eta4, copy='u2', hyper = list(beta = list(fixed = FALSE))) + ",paste(names(covar_z),collapse="+")))



###########

stack_gSG6_test <- inla.stack(
  data = list(y = cbind(as.vector(gsg6cov$pos_gsg6_test), NA, NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, covariate = covar_z))) 

stack_CSP_test <- inla.stack(
  data = list(y = cbind(NA, as.vector(gsg6cov$pos_csp), NA, NA, NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, b.eta2=1:N)))

stack_u1_test <- inla.stack(
  data = list(y = cbind(NA, NA, rep(0,N), NA, NA)),
  A = list(A,diag(N)),
  effects = list(list(field_gSG6 = 1:spde$n.spde),list(intercept_gSG6 = 1, u1 = 1:N, w1 = -1, covariate = covar_z)))

stack_PCR_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, as.vector(gsg6cov$pos_pcr), NA)),
  A = list(A,diag(N)), 
  effects = list(list(field_PCR = 1:spde$n.spde),list(intercept_PCR = 1, b.eta3=1:N, b.eta4=1:N)))

stack_u2_test <- inla.stack(
  data = list(y = cbind(NA, NA, NA, NA, rep(0,N))),
  A = list(A,diag(N)),
  effects = list(list(field_CSP = 1:spde$n.spde),list(intercept_CSP = 1, u2 = 1:N, w2 = -1)))


stack_test <- inla.stack(stack_gSG6_test, stack_CSP_test, stack_u1_test, stack_PCR_test, stack_u2_test)


joint_gsg6_tt20 <- inla(form, c("binomial","binomial","gaussian","binomial","gaussian"),
                        Ntrial = gsg6cov$Nscreen, 
                        data = inla.stack.data(stack_test), 
                        control.predictor = list(A = inla.stack.A(stack_test),compute=TRUE),
                        control.inla = list(int.strategy = 'eb'),
                        control.compute = list(config = TRUE),
                        control.family = list(list(), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE))), list(), list(hyper = list(prec = list(initial = 10, fixed=TRUE)))))


summary(joint_gsg6_tt20)


pred_gSG6_mean_tt20 <- joint_gsg6_tt20$summary.fitted.values[1:N,1]
pred_CSP_mean_tt20 <- joint_gsg6_tt20$summary.fitted.values[(N+1):(2*N),1]

pred_PCR_mean_tt20 <- joint_gsg6_tt20$summary.fitted.values[((3*N)+1):(4*N),1]

pred_gSG6_sd_tt20 <- joint_gsg6_tt20$summary.fitted.values[1:N,2]
pred_CSP_sd_tt20 <- joint_gsg6_tt20$summary.fitted.values[(N+1):(2*N),2]
pred_PCR_sd_tt20 <- joint_gsg6_tt20$summary.fitted.values[((3*N)+1):(4*N),2]

pred_gSG6_ll_tt20 <- joint_gsg6_tt20$summary.fitted.values[1:N,3]
pred_CSP_ll_tt20 <- joint_gsg6_tt20$summary.fitted.values[(N+1):(2*N),3]
pred_PCR_ll_tt20 <- joint_gsg6_tt20$summary.fitted.values[((3*N)+1):(4*N),3]

pred_gSG6_ul_tt20 <- joint_gsg6_tt20$summary.fitted.values[1:N,5]
pred_CSP_ul_tt20 <- joint_gsg6_tt20$summary.fitted.values[(N+1):(2*N),5]
pred_PCR_ul_tt20 <- joint_gsg6_tt20$summary.fitted.values[((3*N)+1):(4*N),5]

pred_tt20 <- cbind(pred_gSG6_mean_tt20, pred_gSG6_sd_tt20, pred_gSG6_ll_tt20, pred_gSG6_ul_tt20, 
                   pred_CSP_mean_tt20 , pred_CSP_sd_tt20, pred_CSP_ll_tt20 , pred_CSP_ul_tt20, 
                   pred_PCR_mean_tt20, pred_PCR_sd_tt20,  pred_PCR_ll_tt20, pred_PCR_ul_tt20)
pred_tt20 <-pred_tt20[test_ind,]
obs_tt20 <- gsg6cov[test_ind,c(1:4,7)]


obs_tt <- rbind(obs_tt1,obs_tt2,obs_tt3,obs_tt4,obs_tt5,obs_tt6,obs_tt7,obs_tt8,obs_tt9,obs_tt10,
                obs_tt11,obs_tt12,obs_tt13,obs_tt14,obs_tt15,obs_tt16,obs_tt17,obs_tt18,obs_tt19,obs_tt20)
pred_tt <- rbind(pred_tt1,pred_tt2,pred_tt3,pred_tt4,pred_tt5,pred_tt6,pred_tt7,pred_tt8,pred_tt9,pred_tt10,
                 pred_tt11,pred_tt12,pred_tt13,pred_tt14,pred_tt15,pred_tt16,pred_tt17,pred_tt18,pred_tt19,pred_tt20)



pred_dat <- data.frame(obs_tt,
                       pred_gsg6  = pred_tt[,1], pred_gsg6_sd=pred_tt[,2], pred_gsg6_ll=pred_tt[,3], pred_gsg6_ul=pred_tt[,4],
                       pred_csp  = pred_tt[,5], pred_csp_sd=pred_tt[,6], pred_csp_ll=pred_tt[,7], pred_csp_ul=pred_tt[,8], 
                       pred_pcr  = pred_tt[,9], pred_pcr_sd=pred_tt[,10], pred_pcr_ll=pred_tt[,11], pred_pcr_ul=pred_tt[,12])

#calculate obs med, lci + uci prob
gsg6_ll_qbeta <- qbeta(0.025, shape1 = 1+obs_tt$pos_gsg6, shape2 = 1+(obs_tt$Nscreen-obs_tt$pos_gsg6))
gsg6_ul_qbeta <- qbeta(0.975, shape1 = 1+obs_tt$pos_gsg6, shape2 = 1+(obs_tt$Nscreen-obs_tt$pos_gsg6))
gsg6_med_qbeta <- qbeta(0.5, shape1 = 1+obs_tt$pos_gsg6, shape2 = 1+(obs_tt$Nscreen-obs_tt$pos_gsg6))

csp_ll_qbeta <- qbeta(0.025, shape1 = 1+obs_tt$pos_csp, shape2 = 1+(obs_tt$Nscreen-obs_tt$pos_csp))
csp_ul_qbeta <- qbeta(0.975, shape1 = 1+obs_tt$pos_csp, shape2 = 1+(obs_tt$Nscreen-obs_tt$pos_csp))
csp_med_qbeta <- qbeta(0.5, shape1 = 1+obs_tt$pos_csp, shape2 = 1+(obs_tt$Nscreen-obs_tt$pos_csp))

pcr_ll_qbeta <- qbeta(0.025, shape1 = 1+obs_tt$pos_pcr, shape2 = 1+(obs_tt$Nscreen-obs_tt$pos_pcr))
pcr_ul_qbeta <- qbeta(0.975, shape1 = 1+obs_tt$pos_pcr, shape2 = 1+(obs_tt$Nscreen-obs_tt$pos_pcr))
pcr_med_qbeta <- qbeta(0.5, shape1 = 1+obs_tt$pos_pcr, shape2 = 1+(obs_tt$Nscreen-obs_tt$pos_pcr))


full_tt <- data.frame(obs_gsg6  = gsg6_med_qbeta, obs_gsg6_ll=gsg6_ll_qbeta, obs_gsg6_ul=gsg6_ul_qbeta, 
                      obs_csp  = csp_med_qbeta, obs_csp_ll=csp_ll_qbeta, obs_csp_ul=csp_ul_qbeta, 
                      obs_pcr  = pcr_med_qbeta, obs_pcr_ll=pcr_ll_qbeta, obs_pcr_ul=pcr_ul_qbeta, 
                      pred_gsg6  = pred_tt[,1], pred_gsg6_ll=pred_tt[,3], pred_gsg6_ul=pred_tt[,4],
                      pred_csp  = pred_tt[,5], pred_csp_ll=pred_tt[,6], pred_csp_ul=pred_tt[,8], 
                      pred_pcr  = pred_tt[,9], pred_pcr_ll=pred_tt[,11], pred_pcr_ul=pred_tt[,12])

##adjust so no prob<0 or >1
full_tt <- full_tt %>% mutate(pred_gsg6_ll_adj=case_when( pred_gsg6_ll <0 ~ 0,
                                                          pred_gsg6_ll >= 0 ~ pred_gsg6_ll))


full_tt <- full_tt %>% mutate(pred_gsg6_adj=case_when( pred_gsg6 >1 ~ 1,
                                                       pred_gsg6 <0 ~ 0,
                                                       pred_gsg6 <= 1 ~ pred_gsg6))

full_tt <- full_tt %>% mutate(pred_gsg6_ul_adj=case_when( pred_gsg6_ul >1 ~ 1,
                                                          pred_gsg6_ul <0 ~ 0,
                                                          pred_gsg6_ul <= 1 ~ pred_gsg6_ul))

full_tt <- full_tt %>% mutate(pred_csp_ll_adj=case_when( pred_csp_ll <0 ~ 0,
                                                           pred_csp_ll >= 0 ~ pred_csp_ll))


full_tt <- full_tt %>% mutate(pred_csp_adj=case_when( pred_csp >1 ~ 1,
                                                        pred_csp <0 ~ 0,
                                                        pred_csp <= 1 ~ pred_csp))

full_tt <- full_tt %>% mutate(pred_csp_ul_adj=case_when( pred_csp_ul >1 ~ 1,
                                                           pred_csp_ul <0 ~ 0,
                                                           pred_csp_ul <= 1 ~ pred_csp_ul))

full_tt <- full_tt %>% mutate(pred_pcr_ll_adj=case_when( pred_pcr_ll <0 ~ 0,
                                                           pred_pcr_ll >= 0 ~ pred_pcr_ll))


full_tt <- full_tt %>% mutate(pred_pcr_adj=case_when( pred_pcr >1 ~ 1,
                                                        pred_pcr <0 ~ 0,
                                                        pred_pcr <= 1 ~ pred_pcr))

full_tt <- full_tt %>% mutate(pred_pcr_ul_adj=case_when( pred_pcr_ul >1 ~ 1,
                                                           pred_pcr_ul <0 ~ 0,
                                                           pred_pcr_ul <= 1 ~ pred_pcr_ul))


pred_dat <- pred_dat %>% mutate(pred_gsg6_adj=case_when( pred_gsg6 >1 ~ 1,
                                                         pred_gsg6 <0 ~ 0,
                                                         pred_gsg6 <= 1 ~ pred_gsg6))

pred_dat <- pred_dat %>% mutate(pred_csp_adj=case_when( pred_csp >1 ~ 1,
                                                          pred_csp <0 ~ 0,
                                                          pred_csp <= 1 ~ pred_csp))

pred_dat <- pred_dat %>% mutate(pred_pcr_adj=case_when( pred_pcr >1 ~ 1,
                                                          pred_pcr <0 ~ 0,
                                                          pred_pcr <= 1 ~ pred_pcr))


#####
pred_dat <- pred_dat %>% mutate(pr_csp =pos_csp/Nscreen)
pred_dat <- pred_dat %>% mutate(pr_pcr =pos_pcr/Nscreen)
###########################################
#gSG6 pred v Obs - binned on pred
#################################
pred_dat <- pred_dat %>% mutate(gsg6_bin =ntile(pred_gsg6_adj, n=10))
pred_dat <- pred_dat %>% mutate(k1 =1)
pred_dat <- within(pred_dat, gsg6bin_no_villages <- ave(gsg6_bin, list(gsg6_bin, k1), FUN=length))

head(pred_dat)

pred_dat <- pred_dat %>% 
  group_by(gsg6_bin) %>% 
  mutate(gsg6bin_sum_n = sum(Nscreen)) %>% 
  ungroup()

pred_dat <- pred_dat %>% 
  group_by(gsg6_bin) %>% 
  mutate(gsg6bin_sum_gsg6 = sum(pos_gsg6)) %>% 
  ungroup()


#################################
#gsg6 percentile weights
pred_dat <- pred_dat %>% mutate(pred_dat,mean_gsg6_sg6pctweight=gsg6bin_no_villages*pred_gsg6_adj*(Nscreen/gsg6bin_sum_n))
pred_dat <- pred_dat %>% mutate(pred_dat,sd_sq_gsg6_sg6pctweight=gsg6bin_no_villages*pred_gsg6_sd^2*(Nscreen/gsg6bin_sum_n))


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

agg_tt <- data.frame(obs_gsg6=gsg6_med_qbeta_agg, obs_gsg6_ll=gsg6_ll_qbeta_agg,obs_gsg6_ul=gsg6_ul_qbeta_agg,
                     pred_gsg6=pred_sg6pct_weights$mean_gsg6_sg6pctweight, pred_gsg6_ll=pred_sg6pct_weights$ll_gsg6_sg6pctweight,pred_gsg6_ul=pred_sg6pct_weights$ul_gsg6_sg6pctweight)

##adjust so no prob<0 or >1
agg_tt <- agg_tt %>% mutate(pred_gsg6_ll_adj=case_when( pred_gsg6_ll <0 ~ 0,
                                                        pred_gsg6_ll >= 0 ~ pred_gsg6_ll))


agg_tt <- agg_tt %>% mutate(pred_gsg6_adj=case_when( pred_gsg6 >1 ~ 1,
                                                     pred_gsg6 <0 ~ 0,
                                                     pred_gsg6 <= 1 ~ pred_gsg6))

agg_tt <- agg_tt %>% mutate(pred_gsg6_ul_adj=case_when( pred_gsg6_ul >1 ~ 1,
                                                        pred_gsg6_ul <0 ~ 0,
                                                        pred_gsg6_ul <= 1 ~ pred_gsg6_ul))



ggplot() +
  geom_point(data = full_tt, mapping =aes(x = pred_gsg6_adj, y = obs_gsg6), col="pink") +
  geom_pointrange(data = full_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6, xmin = pred_gsg6_ll, xmax = pred_gsg6_ul), col="pink")+
  geom_pointrange(data = full_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6, ymin = obs_gsg6_ll, ymax = obs_gsg6_ul), col="pink") +
  geom_point(data = agg_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6), alpha=0.8, col="black") +
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6, xmin = pred_gsg6_ll_adj, xmax = pred_gsg6_ul_adj), alpha=0.8, col="black")+
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6, ymin = obs_gsg6_ll, ymax = obs_gsg6_ul), alpha=0.8, col="black") +
  scale_y_continuous("Observed", breaks = seq(from = 0, to = 1, by = .1), limits=c(0,1)) +
  scale_x_continuous("Predicted", breaks = seq(from = 0,to = 1,by = .10), limits=c(0,1)) +
  theme_classic()+ggtitle("gSG6 spatial model fit")  +
  geom_smooth(data = agg_tt, mapping =aes(x = pred_gsg6_adj, y = obs_gsg6), col="black", linetype = "dashed", method = "lm", se=FALSE, fullrange=TRUE)

cor(agg_tt$obs_gsg6,agg_tt$pred_gsg6_adj)
cor(full_tt$obs_gsg6,full_tt$pred_gsg6)

ggsave(
  filename = "H://Documents/test data/03 Output/gsg6_csp_pcr/predVobs_gsg6.png",
  plot = last_plot(),
  device = NULL,
  scale = 1,
  width = 20,
  height = 15,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

###########################################
#CSP pred v Obs - binned on pred
#################################
pred_dat <- pred_dat %>% mutate(csp_bin =ntile(pred_csp_adj, n=10))
pred_dat <- pred_dat %>% mutate(k1 =1)
pred_dat <- within(pred_dat, cspbin_no_villages <- ave(csp_bin, list(csp_bin, k1), FUN=length))

head(pred_dat)

pred_dat <- pred_dat %>% 
  group_by(csp_bin) %>% 
  mutate(cspbin_sum_n = sum(Nscreen)) %>% 
  ungroup()

pred_dat <- pred_dat %>% 
  group_by(csp_bin) %>% 
  mutate(cspbin_sum_csp = sum(pos_csp)) %>% 
  ungroup()


#################################
#csp percentile weights
pred_dat <- pred_dat %>% mutate(pred_dat,mean_csp_sg6pctweight=cspbin_no_villages*pred_csp_adj*(Nscreen/cspbin_sum_n))
pred_dat <- pred_dat %>% mutate(pred_dat,sd_sq_csp_sg6pctweight=cspbin_no_villages*pred_csp_sd^2*(Nscreen/cspbin_sum_n))


pred_sg6pct_weights = aggregate(pred_dat,
                                by = list(pred_dat$csp_bin),
                                FUN = mean)

pred_sg6pct_weights$sd_csp_sg6pctweight  <- sqrt(pred_sg6pct_weights$sd_sq_csp_sg6pctweight)


head(pred_sg6pct_weights)
##calculate 95%CredI

pred_sg6pct_weights <- pred_sg6pct_weights %>% mutate(pred_sg6pct_weights,ll_csp_sg6pctweight=mean_csp_sg6pctweight-(2*sd_csp_sg6pctweight))
pred_sg6pct_weights <- pred_sg6pct_weights %>% mutate(pred_sg6pct_weights,ul_csp_sg6pctweight=mean_csp_sg6pctweight+(2*sd_csp_sg6pctweight))


csp_ll_qbeta_agg <- qbeta(0.025, shape1 = 1+pred_sg6pct_weights$cspbin_sum_csp, shape2 = 1+(pred_sg6pct_weights$cspbin_sum_n-pred_sg6pct_weights$cspbin_sum_csp))
csp_ul_qbeta_agg <- qbeta(0.975, shape1 = 1+pred_sg6pct_weights$cspbin_sum_csp, shape2 = 1+(pred_sg6pct_weights$cspbin_sum_n-pred_sg6pct_weights$cspbin_sum_csp))
csp_med_qbeta_agg <- qbeta(0.5, shape1 = 1+pred_sg6pct_weights$cspbin_sum_csp, shape2 = 1+(pred_sg6pct_weights$cspbin_sum_n-pred_sg6pct_weights$cspbin_sum_csp))

agg_tt <- data.frame(agg_tt, obs_csp=csp_med_qbeta_agg, obs_csp_ll=csp_ll_qbeta_agg,obs_csp_ul=csp_ul_qbeta_agg,
                     pred_csp=pred_sg6pct_weights$mean_csp_sg6pctweight, pred_csp_ll=pred_sg6pct_weights$ll_csp_sg6pctweight,pred_csp_ul=pred_sg6pct_weights$ul_csp_sg6pctweight)

##adjust so no prob<0 or >1
agg_tt <- agg_tt %>% mutate(pred_csp_ll_adj=case_when( pred_csp_ll <0 ~ 0,
                                                         pred_csp_ll >= 0 ~ pred_csp_ll))


agg_tt <- agg_tt %>% mutate(pred_csp_adj=case_when( pred_csp >1 ~ 1,
                                                      pred_csp <0 ~ 0,
                                                      pred_csp <= 1 ~ pred_csp))

agg_tt <- agg_tt %>% mutate(pred_csp_ul_adj=case_when( pred_csp_ul >1 ~ 1,
                                                         pred_csp_ul <0 ~ 0,
                                                         pred_csp_ul <= 1 ~ pred_csp_ul))


ggplot() +
  geom_point(data = full_tt, mapping =aes(x = pred_csp_adj, y = obs_csp), col="pink") +
  geom_pointrange(data = full_tt, mapping = aes(x = pred_csp_adj, y = obs_csp, xmin = pred_csp_ll, xmax = pred_csp_ul), col="pink")+
  geom_pointrange(data = full_tt, mapping = aes(x = pred_csp_adj, y = obs_csp, ymin = obs_csp_ll, ymax = obs_csp_ul), col="pink") +
  geom_point(data = agg_tt, mapping = aes(x = pred_csp_adj, y = obs_csp), alpha=0.8, col="black") +
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_csp_adj, y = obs_csp, xmin = pred_csp_ll_adj, xmax = pred_csp_ul_adj), alpha=0.8, col="black")+
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_csp_adj, y = obs_csp, ymin = obs_csp_ll, ymax = obs_csp_ul), alpha=0.8, col="black") +
  scale_y_continuous("Observed", breaks = seq(from = 0, to = .8, by = .1), limits=c(0,.8)) +
  scale_x_continuous("Predicted", breaks = seq(from = 0,to = .8,by = .10), limits=c(0,.8)) +
  theme_classic()+ggtitle("csp spatial model fit")  +
  geom_smooth(data = agg_tt, mapping =aes(x = pred_csp_adj, y = obs_csp), col="black", linetype = "dashed", method = "lm", se=FALSE, fullrange=TRUE)

cor(agg_tt$obs_csp,agg_tt$pred_csp_adj)
cor(full_tt$obs_csp,full_tt$pred_csp)

ggsave(
  filename = "H://Documents/test data/03 Output/gsg6_csp_pcr/predVobs_csp.png",
  plot = last_plot(),
  device = NULL,
  scale = 1,
  width = 20,
  height = 15,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

###########################################
#PCR pred v Obs - binned on pred
#################################
pred_dat <- pred_dat %>% mutate(pcr_bin =ntile(pred_pcr_adj, n=10))
pred_dat <- pred_dat %>% mutate(k1 =1)
pred_dat <- within(pred_dat, pcrbin_no_villages <- ave(pcr_bin, list(pcr_bin, k1), FUN=length))

head(pred_dat)

pred_dat <- pred_dat %>% 
  group_by(pcr_bin) %>% 
  mutate(pcrbin_sum_n = sum(Nscreen)) %>% 
  ungroup()

pred_dat <- pred_dat %>% 
  group_by(pcr_bin) %>% 
  mutate(pcrbin_sum_pcr = sum(pos_pcr)) %>% 
  ungroup()


#################################
#pcr percentile weights
pred_dat <- pred_dat %>% mutate(pred_dat,mean_pcr_sg6pctweight=pcrbin_no_villages*pred_pcr_adj*(Nscreen/pcrbin_sum_n))
pred_dat <- pred_dat %>% mutate(pred_dat,sd_sq_pcr_sg6pctweight=pcrbin_no_villages*pred_pcr_sd^2*(Nscreen/pcrbin_sum_n))


pred_sg6pct_weights = aggregate(pred_dat,
                                by = list(pred_dat$pcr_bin),
                                FUN = mean)

pred_sg6pct_weights$sd_pcr_sg6pctweight  <- sqrt(pred_sg6pct_weights$sd_sq_pcr_sg6pctweight)


head(pred_sg6pct_weights)
##calculate 95%CredI

pred_sg6pct_weights <- pred_sg6pct_weights %>% mutate(pred_sg6pct_weights,ll_pcr_sg6pctweight=mean_pcr_sg6pctweight-(2*sd_pcr_sg6pctweight))
pred_sg6pct_weights <- pred_sg6pct_weights %>% mutate(pred_sg6pct_weights,ul_pcr_sg6pctweight=mean_pcr_sg6pctweight+(2*sd_pcr_sg6pctweight))


pcr_ll_qbeta_agg <- qbeta(0.025, shape1 = 1+pred_sg6pct_weights$pcrbin_sum_pcr, shape2 = 1+(pred_sg6pct_weights$pcrbin_sum_n-pred_sg6pct_weights$pcrbin_sum_pcr))
pcr_ul_qbeta_agg <- qbeta(0.975, shape1 = 1+pred_sg6pct_weights$pcrbin_sum_pcr, shape2 = 1+(pred_sg6pct_weights$pcrbin_sum_n-pred_sg6pct_weights$pcrbin_sum_pcr))
pcr_med_qbeta_agg <- qbeta(0.5, shape1 = 1+pred_sg6pct_weights$pcrbin_sum_pcr, shape2 = 1+(pred_sg6pct_weights$pcrbin_sum_n-pred_sg6pct_weights$pcrbin_sum_pcr))

agg_tt <- data.frame(agg_tt,obs_pcr=pcr_med_qbeta_agg, obs_pcr_ll=pcr_ll_qbeta_agg,obs_pcr_ul=pcr_ul_qbeta_agg,
                     pred_pcr=pred_sg6pct_weights$mean_pcr_sg6pctweight, pred_pcr_ll=pred_sg6pct_weights$ll_pcr_sg6pctweight,pred_pcr_ul=pred_sg6pct_weights$ul_pcr_sg6pctweight)

##adjust so no prob<0 or >1
agg_tt <- agg_tt %>% mutate(pred_pcr_ll_adj=case_when( pred_pcr_ll <0 ~ 0,
                                                         pred_pcr_ll >= 0 ~ pred_pcr_ll))


agg_tt <- agg_tt %>% mutate(pred_pcr_adj=case_when( pred_pcr >1 ~ 1,
                                                      pred_pcr <0 ~ 0,
                                                      pred_pcr <= 1 ~ pred_pcr))

agg_tt <- agg_tt %>% mutate(pred_pcr_ul_adj=case_when( pred_pcr_ul >1 ~ 1,
                                                         pred_pcr_ul <0 ~ 0,
                                                         pred_pcr_ul <= 1 ~ pred_pcr_ul))



ggplot() +
  geom_point(data = full_tt, mapping =aes(x = pred_pcr_adj, y = obs_pcr), col="pink") +
  geom_pointrange(data = full_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr, xmin = pred_pcr_ll, xmax = pred_pcr_ul), col="pink")+
  geom_pointrange(data = full_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr, ymin = obs_pcr_ll, ymax = obs_pcr_ul), col="pink") +
  geom_point(data = agg_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr), alpha=0.8, col="black") +
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr, xmin = pred_pcr_ll_adj, xmax = pred_pcr_ul_adj), alpha=0.8, col="black")+
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr, ymin = obs_pcr_ll, ymax = obs_pcr_ul), alpha=0.8, col="black") +
  scale_y_continuous("Observed", breaks = seq(from = 0, to = .3, by = .05), limits=c(0,.3)) +
  scale_x_continuous("Predicted", breaks = seq(from = 0,to = .1,by = .05), limits=c(0,.1)) +
  theme_classic()+ggtitle("pcr spatial model fit")  +
  geom_smooth(data = agg_tt, mapping =aes(x = pred_pcr_adj, y = obs_pcr), col="black", linetype = "dashed", method = "lm", se=FALSE, fullrange=TRUE)

cor(agg_tt$obs_pcr,agg_tt$pred_pcr_adj)
cor(full_tt$obs_pcr,full_tt$pred_pcr)

ggsave(
  filename = "H://Documents/test data/03 Output/gsg6_csp_pcr/predVobs_pcr.png",
  plot = last_plot(),
  device = NULL,
  scale = 1,
  width = 20,
  height = 15,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

##################################
agg_tt <- agg_tt*100
full_tt <- full_tt*100


ggplot() +
  geom_point(data = full_tt, mapping =aes(x = pred_gsg6_adj, y = obs_gsg6), col="pink") +
  geom_pointrange(data = full_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6, xmin = pred_gsg6_ll, xmax = pred_gsg6_ul), col="pink")+
  geom_pointrange(data = full_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6, ymin = obs_gsg6_ll, ymax = obs_gsg6_ul), col="pink") +
  geom_point(data = agg_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6), alpha=0.8, col="black") +
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6, xmin = pred_gsg6_ll_adj, xmax = pred_gsg6_ul_adj), alpha=0.8, col="black")+
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_gsg6_adj, y = obs_gsg6, ymin = obs_gsg6_ll, ymax = obs_gsg6_ul), alpha=0.8, col="black") +
  scale_y_continuous("Observed gSG6 Prevalence", breaks = seq(from = 0, to = 100, by = 20), limits=c(0,100)) +
  scale_x_continuous("Predicted gSG6 Prevalence", breaks = seq(from = 0,to = 100,by = 20), limits=c(0,100)) +
  theme_classic()+ggtitle("gSG6 spatial model fit")  +
  geom_smooth(data = agg_tt, mapping =aes(x = pred_gsg6_adj, y = obs_gsg6), col="black", linetype = "dashed", method = "lm", se=FALSE, fullrange=TRUE)

cor(agg_tt$obs_gsg6,agg_tt$pred_gsg6_adj)
cor(full_tt$obs_gsg6,full_tt$pred_gsg6)

ggsave(
  filename = "H://Documents/test data/03 Output/gsg6_csp_pcr/predVobs_gsg6_prev.png",
  plot = last_plot(),
  device = NULL,
  scale = 1,
  width = 20,
  height = 15,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)


ggplot() +
  geom_point(data = full_tt, mapping =aes(x = pred_csp_adj, y = obs_csp), col="pink") +
  geom_pointrange(data = full_tt, mapping = aes(x = pred_csp_adj, y = obs_csp, xmin = pred_csp_ll, xmax = pred_csp_ul), col="pink")+
  geom_pointrange(data = full_tt, mapping = aes(x = pred_csp_adj, y = obs_csp, ymin = obs_csp_ll, ymax = obs_csp_ul), col="pink") +
  geom_point(data = agg_tt, mapping = aes(x = pred_csp_adj, y = obs_csp), alpha=0.8, col="black") +
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_csp_adj, y = obs_csp, xmin = pred_csp_ll_adj, xmax = pred_csp_ul_adj), alpha=0.8, col="black")+
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_csp_adj, y = obs_csp, ymin = obs_csp_ll, ymax = obs_csp_ul), alpha=0.8, col="black") +
  scale_y_continuous("Observed csp Prevalence", breaks = seq(from = 0, to = 80, by = 20), limits=c(0,80)) +
  scale_x_continuous("Predicted csp Prevalence", breaks = seq(from = 0,to = 80,by = 20), limits=c(0,80)) +
  theme_classic()+ggtitle("csp spatial model fit")  +
  geom_smooth(data = agg_tt, mapping =aes(x = pred_csp_adj, y = obs_csp), col="black", linetype = "dashed", method = "lm", se=FALSE, fullrange=TRUE)

cor(agg_tt$obs_csp,agg_tt$pred_csp_adj)
cor(full_tt$obs_csp,full_tt$pred_csp)

ggsave(
  filename = "H://Documents/test data/03 Output/gsg6_csp_pcr/predVobs_csp_prev.png",
  plot = last_plot(),
  device = NULL,
  scale = 1,
  width = 20,
  height = 15,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)



ggplot() +
  geom_point(data = full_tt, mapping =aes(x = pred_pcr_adj, y = obs_pcr), col="pink") +
  geom_pointrange(data = full_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr, xmin = pred_pcr_ll, xmax = pred_pcr_ul), col="pink")+
  geom_pointrange(data = full_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr, ymin = obs_pcr_ll, ymax = obs_pcr_ul), col="pink") +
  geom_point(data = agg_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr), alpha=0.8, col="black") +
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr, xmin = pred_pcr_ll_adj, xmax = pred_pcr_ul_adj), alpha=0.8, col="black")+
  geom_pointrange(data = agg_tt, mapping = aes(x = pred_pcr_adj, y = obs_pcr, ymin = obs_pcr_ll, ymax = obs_pcr_ul), alpha=0.8, col="black") +
  scale_y_continuous("Observed pcr Prevalence", breaks = seq(from = 0, to = 30, by = 5), limits=c(0,30)) +
  scale_x_continuous("Predicted pcr Prevalence", breaks = seq(from = 0,to = 8,by = 2), limits=c(0,8)) +
  theme_classic()+ggtitle("pcr spatial model fit")  +
  geom_smooth(data = agg_tt, mapping =aes(x = pred_pcr_adj, y = obs_pcr), col="black", linetype = "dashed", method = "lm", se=FALSE, fullrange=TRUE)

cor(agg_tt$obs_pcr,agg_tt$pred_pcr_adj)
cor(full_tt$obs_pcr,full_tt$pred_pcr)

ggsave(
  filename = "H://Documents/test data/03 Output/gsg6_csp_pcr/predVobs_pcr_prev.png",
  plot = last_plot(),
  device = NULL,
  scale = 1,
  width = 20,
  height = 15,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)
