setwd("~/RT-Shrinkage-Paper-Codes/data sets")
rm(list=ls()) # remove all files
library(lars) # for 'cv.folds' function to split data
library(ncvreg)
library(glmnet)
library(gdata)
library(readxl)
library(caret)
set.seed(2024) # set seed number for reproducibility
##

df <- readRDS("pollution.rds")
Xr <- as.matrix(df_file$X)
yr <- as.matrix(df_file$y)
##


KL_beta <- function(G,k,Ip,beta_hat){
  return(solve(G+k*Ip)%*%(G-k*Ip)%*%beta_hat)
}

Kl_beta_k_opt <- function(k_grid,X.test,y.test,mean.y,mean.X,sd.X,beta_hat){
  MSE_beta <- c()
  p <- dim(X.test)[2]
  for (i in 1:length(k_grid)) {
    beta.KL <- KL_beta(G,k_grid[i],Ip,beta_hat)
    MSE_beta[i] <- PE.funct_beta(beta.KL,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  }
  k_opt = k_grid[which.min(MSE_beta)]
  beta_opt = KL_beta(G,k_opt,Ip,beta_hat)
  MSE = MSE_beta[which.min(MSE_beta)]
  #plot(MSE_beta)
  return(list(k_opt = k_opt,
              beta_opt = beta_opt,
              MSE = MSE))
}

LS_beta_d_opt <- function(d_grid,X.test,y.test,mean.y,mean.X,sd.X,beta_fm, beta_sm){
  MSE_beta <- c()
  p <- dim(X.test)[2]
  for (i in 1:length(d_grid)) {
    beta.LS <- beta_fm - d_grid[i]*(beta_fm-beta_sm)
    MSE_beta[i] <- PE.funct_beta(beta.LS,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  }
  d_opt = d_grid[which.min(MSE_beta)]
  beta_opt = beta_fm - d_opt*(beta_fm-beta_sm)
  MSE = MSE_beta[which.min(MSE_beta)]
  #plot(MSE_beta)
  return(list(d_opt = d_opt,
              beta_opt = beta_opt,
              MSE = MSE))
}

SPTKL_beta_d_opt <- function(d_grid,X.test,y.test,mean.y,mean.X,sd.X,beta_fm, beta_sm,tn,lalpha){
  MSE_beta <- c()
  p <- dim(X.test)[2]
  for (i in 1:length(d_grid)) {
    beta.SPTKL  <- beta_fm - d_grid[i]*(beta_fm-beta_sm)*I(tn < lalpha)
    MSE_beta[i] <- PE.funct_beta(beta.SPTKL,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  }
  d_opt = d_grid[which.min(MSE_beta)]
  beta_opt = beta_fm - d_opt*(beta_fm-beta_sm)*I(tn < lalpha)
  MSE = MSE_beta[which.min(MSE_beta)]
  #plot(MSE_beta)
  return(list(d_opt = d_opt,
              beta_opt = beta_opt,
              MSE = MSE))
}

################################################################

################################################################
Data.matrix <- data.frame(yr,Xr)
# n and p parameters
n <- dim(Xr)[1]
p <- dim(Xr)[2]
#
full_ind <- colnames(Xr)
# Variable Selection
require(olsrr)
model <- lm(yr ~ ., data = Data.matrix)
bss <- ols_step_best_subset(model)
#### Sub Model Variables based onn CP and AIC Criteria
sub1_x <- unlist(strsplit(bss$metrics$predictors[which.min(bss$metrics$cp)], " "))
sub2_x <- unlist(strsplit(bss$metrics$predictors[which.min(bss$metrics$aic)], " "))
#
p1_indx_1 <- match(sub1_x,colnames(Xr),nomatch = NA_integer_)
p1_indx_2 <- match(sub2_x,colnames(Xr),nomatch = NA_integer_)
### start here
p1_1 <- length(p1_indx_1)
p1_2 <- length(p1_indx_2)

p2_1 <- p-p1_1
p2_2 <- p-p1_2
## Unrestricted Model Formulation without the intercept term
xcount.ur <- c(0,full_ind) 
fmla.ur <- as.formula(paste("y.train ~ ", paste(xcount.ur, collapse= "+")))

# K-fold cross-validation
### PREDICTION ERROR (You may use different function if you want)
# This is constucted with standardized train data 
PE.funct_beta <- function(opt.beta, newx, p1_indx,newy,mean.y,mean.X,sd.X)
{
  fit = newy - mean.y -  scale(newx[,p1_indx], mean.X[p1_indx],F) %*% opt.beta
  return(t(fit)%*%fit / length(newy))
}

#PE.funct_beta <- function(opt.beta, newx, p1_indx,newy,mean.y,mean.X)
#{
#  fit = newy - newx[,p1_indx] %*% opt.beta
#  return(t(fit)%*%fit / length(newy))
#}

prederror<-function(y,yhat){
  mean((yhat - y)^2)
}

#######
# Create a resampling method
trainControl_all <- trainControl(
  method = "cv", 
  number = 5
)
## network grids tuning
nnetGrid <-  expand.grid(size = seq(from = 1, to = 10, by = 1),
                         decay = seq(from = 0.1, to = 0.5, by = 0.1))

####### Whether the TRAIN data is standardized or not
standardize = TRUE

# The number of bootstrap replication

repeatnum <- 1000



#### SAVE MATRICEs, you may increase the number of ncol along with the numbers of estimators
coltitles.MSE <- c('KL','RKL','KLPT','KLS','KLPS','LSKL','SPTKL')

MSE_sub1 <- MSE_sub2 <-MSE_sub3 <-MSE_sub4 <- matrix(NA, nrow=repeatnum,ncol=7)
colnames(MSE_sub1) <-colnames(MSE_sub2) <-colnames(MSE_sub3)  <-colnames(MSE_sub4) <- coltitles.MSE

coltitles.MSE_pens <- c('ENET','LASSO','RIDGE','ALASSO','SCAD')
MSE_pens    <- matrix(NA, nrow=repeatnum,ncol=5)
colnames(MSE_pens) <- coltitles.MSE_pens

coltitles.MSE_ml <- c('NN','RF','KNN','LSE')
MSE_ml    <- matrix(NA, nrow=repeatnum,ncol=4)
colnames(MSE_ml) <- coltitles.MSE_ml

# Grids
Ip=diag(p)
k_grid <- c(10^seq(-1,1,length=25),10^seq(1,2,length=25)[-1])
d_grid <- seq(0,1,.1)


smp_size <- floor(0.7 * n)
#set.seed(2024) # set seed number for reproducibility
for(J in 1:repeatnum){ # loop 1 begins
  set.seed(J)
  #generare data with resampling
  #data.boots <- as.matrix(Data.matrix[sample(1:nrow(Data.matrix), nrow(Data.matrix), replace = TRUE),] )
  #sample <- sample(c(TRUE, FALSE), nrow(prostate), replace=TRUE, prob=c(0.8,0.2))
  train_ind <- sample(seq_len(nrow(Data.matrix)), size = smp_size)
  #prostate_train <- prostate[sample, ]
  #prostate_test <- prostate[!sample, ]
  X.train <- Data.matrix[train_ind,-1]
  y.train <- Data.matrix[train_ind,1]
  nt <- dim(X.train)[1]
  ############################################################################################################################################
  ############################################################################################################################################
  ############################################################################################################################################
 
  
  if (standardize == TRUE) {
    mean.X  <- apply (X.train, 2, mean)
    sd.X    <- apply (X.train, 2, sd)
    X.train <- scale(X.train, mean.X, F)
    mean.y  <- mean (y.train)
    y.train <- y.train - mean.y 
  } else {
    mean.X <- rep (0, length = ncol(X.train)) 
    mean.y <- 0}
  
  
  train_data <- data.frame(y.train,X.train)
  
  X.test  <- Data.matrix[-train_ind,-1]
  y.test  <- Data.matrix[-train_ind,1]
  y.test_scale <- y.test - mean.y
  X.test_scale <- scale(X.test, mean.X, F)
  
  
  
  
  
  # Unrestricted fit
  ufit<-lm(fmla.ur, data =  train_data)
  beta.FM<-ufit$coef
  beta1.FM_1<-beta1.FM_2<-beta1.FM_3<-beta1.FM_4<-beta.FM
  #
  residuals_ufit    = residuals(ufit)
  sigma2_1          = sum(residuals_ufit^2) / (length(residuals_ufit) - p)
  sigma2_2          = sum(residuals_ufit^2) / (length(residuals_ufit) - p)
  
  
  R1 = diag(1,p,p)[-c(p1_indx_1),]
  R2 = diag(1,p,p)[-c(p1_indx_2),]
  G <- t(X.train) %*% (X.train)
  
  
  beta1.SM_1=beta.FM-solve(G)%*%t(R1)%*%solve(R1%*%solve(G)%*%t(R1))%*%R1%*%beta.FM
  beta1.SM_2=beta.FM-solve(G)%*%t(R2)%*%solve(R2%*%solve(G)%*%t(R2))%*%R2%*%beta.FM
  
  
  
  
  
  #sigma(rfit)^2 == sigma2 # Check they are equal or not!
  
  
  
  
  #beta.KL     = solve(G+k*Ip)%*%(G-k*Ip)%*%beta.FM
  beta.KL      = Kl_beta_k_opt(k_grid,X.test,y.test,mean.y,mean.X,sd.X,beta.FM)$beta_opt
  #
  beta.RKL1    = Kl_beta_k_opt(k_grid,X.test,y.test,mean.y,mean.X,sd.X,beta1.SM_1)$beta_opt
  beta.RKL2    = Kl_beta_k_opt(k_grid,X.test,y.test,mean.y,mean.X,sd.X,beta1.SM_2)$beta_opt

  
  
  tn1=drop(t(R1%*%beta.FM)%*%solve(R1%*%solve(G)%*%t(R1))%*%(R1%*%beta.FM)/(p2_1*sigma2_1))
  tn2=drop(t(R2%*%beta.FM)%*%solve(R2%*%solve(G)%*%t(R2))%*%(R2%*%beta.FM)/(p2_2*sigma2_2))
  
  
  piconstan1=(p2_1-2)*(nt-p)/(p2_1*(nt-p+2))
  piconstan2=(p2_2-2)*(nt-p)/(p2_2*(nt-p+2))
  
  
  lalpha1=qf(.95,p2_1,nt-p)
  lalpha2=qf(.95,p2_2,nt-p)
  
  
  beta1.PT_1 <-beta.KL-(beta.KL-beta.RKL1)*I(tn1 < lalpha1 )
  beta1.PT_2 <-beta.KL-(beta.KL-beta.RKL2)*I(tn2 < lalpha2 )

  
  beta1.S_1  <-beta.KL-((beta.KL-beta.RKL1)*piconstan1/tn1)
  beta1.S_2  <-beta.KL-((beta.KL-beta.RKL2)*piconstan2/tn2)

  beta1.PS_1   <-beta.RKL1+ max(0,(1-piconstan1/tn1))*(beta.KL-beta.RKL1)
  beta1.PS_2   <-beta.RKL2+ max(0,(1-piconstan2/tn2))*(beta.KL-beta.RKL2)

  
  beta.LSKL_1   <- LS_beta_d_opt(d_grid,X.test,y.test,mean.y,mean.X,sd.X,beta.KL,beta1.SM_1)$beta_opt
  beta.LSKL_2   <- LS_beta_d_opt(d_grid,X.test,y.test,mean.y,mean.X,sd.X,beta.KL,beta1.SM_2)$beta_opt

  
  
  beta.SPTKL_1   <- SPTKL_beta_d_opt(d_grid,X.test,y.test,mean.y,mean.X,sd.X,beta.KL,beta1.SM_1,tn1,lalpha1)$beta_opt
  beta.SPTKL_2   <- SPTKL_beta_d_opt(d_grid,X.test,y.test,mean.y,mean.X,sd.X,beta.KL,beta1.SM_2,tn2,lalpha2)$beta_opt

  #PE.funct_beta(beta.FM,X.test,1:p,y.test,mean.y,mean.X)
  #prederror(y.test-mean.y,scale (X.test, mean.X, FALSE) %*% beta.FM)
  #prederror(y.test-mean.y,X.test_scale %*% beta.FM)
  #cbind(y.test-mean.y,scale (X.test, mean.X, FALSE) %*% beta.FM, predict(ufit, newdata=as.data.frame(X.test_scale)))
  
  
  #### PENALIZED METHODS
  ##### ENET
  alphas <- seq(.1,.9,.1)
  fits.enet <- list()
  for (ind in 1:length(alphas)){
    fits.enet[[ind]] <- cv.glmnet(X.train, y.train, alpha = alphas[ind],standardize=T)
  }
  cvs.enet <- sapply(fits.enet, function(x) {
    beta.enet = coef(x,s="lambda.min")[-1]
    return(list(PE=mean((y.test_scale - X.test_scale%*%beta.enet)^2),beta.enet = beta.enet))})
  ind_opt  <- which.min(do.call(rbind,cvs.enet[1,]))
  b.enet <- do.call(rbind,cvs.enet[2,])[ind_opt,]
  #  
  lasso.fit <- cv.glmnet(X.train, y.train, alpha = 1,standardize=T)
  b.lasso = coef(lasso.fit,s="lambda.min")[-1]
  
  
  
  
  
  ridge.fit <- cv.glmnet(X.train, y.train, alpha = 0,standardize=T)
  b.ridge = coef(ridge.fit,s="lambda.min")[-1]
  
  
  
  weight <- b.lasso
  weight <- ifelse(weight == 0, 0.00001, weight)
  #adalaso.cv.model <- cv.glmnet(X.train, y.train, alpha = 1, penalty.factor=1/abs(weight),nfolds = 10,intercept = F, standardize = F)
  #beta.alasso   <- coef(adalaso.cv.model, s = "lambda.min")[-1]
  alasso.fit <- cv.glmnet(X.train, y.train, alpha = 1, penalty.factor=1/abs(weight),standardize = T)
  b.alasso = coef(alasso.fit,s="lambda.min")[-1]
  
  
  
  #scad.cv.model <- cv.ncvreg(X.train, y.train, penalty=c("SCAD"),nfolds = 10,intercept = T, standardize = F)
  #b.scad        <- coef(scad.cv.model, s = "lambda.min")[-1]
  
  #SCAD
  scad.fit <- cv.ncvreg(X.train, y.train, penalty=c("SCAD"))
  b.scad   <- coef(scad.fit,s="lambda.min")[-1]
  
  
  
  # Machine Learnings
  
  ##########################################################
  ###### Neurol Network
  ##########################################################
  train_df <- data.frame(y=y.train,X.train)
  test_df  <- data.frame(y=y.test_scale,X.test_scale)
  df_center <- rbind(train_df,test_df)
  train_indx <- 1:dim(train_data)[1]
  #
  max <-  apply(df_center, 2 ,max)
  min <-  apply(df_center, 2 ,min)
  scale_df <-  as.data.frame(scale(df_center, center = min, scale = max - min))
  #
  X.train_nnet   <- scale_df[train_indx,-1]
  y.train_nnet   <- scale_df[train_indx, 1]
  train_data_nnet <- data.frame(y.train_nnet,X.train_nnet)
  #
  X.test_nnet  <- scale_df[-train_indx,-1]
  y.test_nnet  <- scale_df[-train_indx,1]
  
  
  annMod <- train(form = y.train_nnet ~., # use all other variables to predict target
                  data = train_data_nnet, # training data
                  #preProcess = "range", # apply min-max normalization
                  method = "nnet", # use nnet()
                  trControl = trainControl_all,
                  metric = "RMSE",
                  tuneGrid = nnetGrid, # search over the created grid
                  linout = TRUE,
                  trace = FALSE) # suppress output
  
  
  
  ##########################################################
  ##########################################################
  ##      Random Forest
  ##########################################################
  ##########################################################
  m.randomForest <- train(y.train ~ ., 
                          data = train_data, 
                          method = "rf", 
                          trControl = trainControl_all,
                          na.action = na.omit,
                          trace = FALSE)
  
  
  ##########################################################
  ##########################################################
  ##      KNN regression
  ##########################################################
  ##########################################################
  # Create a hyperparameter grid search
  hyper_grid <- expand.grid(
    k = floor(seq(1, nrow(X.train)/3, length.out = 20))
  )
  
  # Fit knn model and perform grid search
  knn_grid <- train(
    form = y.train ~., # use all other variables to predict target
    data = train_data, # training data
    method = "knn", 
    trControl = trainControl_all, 
    tuneGrid = hyper_grid,
    metric = "RMSE",
    preProc = c("center", "scale")
  )
  

  
  #prederror(y.test-mean.y,scale(X.test[,p1_indx_1], mean.X[p1_indx_1],FALSE) %*% beta1.FM_1)
  #MSE of based on sub1
  MSE_sub1[J,1] <- PE.funct_beta(beta.KL,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub1[J,2] <- PE.funct_beta(beta.RKL1,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub1[J,3] <- PE.funct_beta(beta1.PT_1,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub1[J,4] <- PE.funct_beta(beta1.S_1,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub1[J,5] <- PE.funct_beta(beta1.PS_1,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub1[J,6] <- PE.funct_beta(beta.LSKL_1,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub1[J,7] <- PE.funct_beta(beta.SPTKL_1,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  
  
  
  #MSE of based on sub2
  MSE_sub2[J,1] <- PE.funct_beta(beta.KL,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub2[J,2] <- PE.funct_beta(beta.RKL2,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub2[J,3] <- PE.funct_beta(beta1.PT_2,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub2[J,4] <- PE.funct_beta(beta1.S_2,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub2[J,5] <- PE.funct_beta(beta1.PS_2,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub2[J,6] <- PE.funct_beta(beta.LSKL_2,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_sub2[J,7] <- PE.funct_beta(beta.SPTKL_2,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  
  
  #PE.funct_beta(b.enet,X.test,1:p,y.test,mean.y,mean.X)
  #PE.funct_beta(b.lasso,X.test,1:p,y.test,mean.y,mean.X)
  #PE.funct_beta(b.alasso,X.test,1:p,y.test,mean.y,mean.X)
  #PE.funct_beta(b.ridge,X.test,1:p,y.test,mean.y,mean.X)
  #PE.funct_beta(b.scad,X.test,1:p,y.test,mean.y,mean.X)
  #PE.funct_beta(beta.FM,X.test,1:p,y.test,mean.y,mean.X)
  
  #MSE of based on sub3
  MSE_pens[J,1] <- PE.funct_beta(b.enet,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_pens[J,2] <- PE.funct_beta(b.lasso,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_pens[J,3] <- PE.funct_beta(b.alasso,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_pens[J,4] <- PE.funct_beta(b.ridge,X.test,1:p,y.test,mean.y,mean.X,sd.X)
  MSE_pens[J,5] <- PE.funct_beta(b.scad,X.test,1:p,y.test,mean.y,mean.X,sd.X)

  
  #MSE ML
  MSE_ml[J,1]   <- prederror(y.test_nnet*(max(df_center$y)-min(df_center$y))+min(df_center$y), 
                             predict(annMod,X.test_nnet)*(max(df_center$y)-min(df_center$y))+min(df_center$y))  
  MSE_ml[J,2]   <- prederror(y.test_scale,predict(m.randomForest,X.test_scale))
  MSE_ml[J,3]   <- prederror(y.test_scale,predict(knn_grid, newdata=X.test_scale))
  MSE_ml[J,4]   <- prederror(y.test_scale,predict(ufit, newdata=as.data.frame(X.test_scale)))
  
  
  
}#loop 1 ends

#save.image("~/Desktop/ALL_MY_STUDIES/Makaleler/marwan al-momani/real_Data_finalized/KL_realdata_pollution_shrinkage_10cv_1000rep.RData")


#See results of Shrinkage
PE.all1     <- apply(MSE_sub1,2,mean)                
PE.all2     <- apply(MSE_sub2,2,mean)                
PE.se.all1 <- apply(MSE_sub1,2,sd)/sqrt(nrow(MSE_sub1))
PE.se.all2 <- apply(MSE_sub2,2,sd)/sqrt(nrow(MSE_sub2))

#See results of Penalized
PE.pens     <- apply(MSE_pens,2,mean)                
PE.se.pens  <- apply(MSE_pens,2,sd)/sqrt(nrow(MSE_pens))

#See results of ML
PE.mls     <- apply(MSE_ml,2,mean)
PE.se.mls  <- apply(MSE_ml,2,sd)/sqrt(nrow(MSE_ml))


T1    <- cbind(PE.all1,PE.all2)
T1_se <- cbind(PE.se.all1,PE.se.all2)

T1[1,1]/T1
T1[1,1]/PE.pens
T1[1,1]/PE.mls
####

T2      <- PE.pens
T2_se   <- PE.se.pens

T1[1,1]/T2

T3      <- PE.mls
T3_se   <- PE.se.mls

T1[1,1]/T3

Table <- round(T1,4)
colnames(T1) <- c('PE','PE')
rownames(T1) <- names(PE.all1)
for(i in 1:7){
  for(j in c(1:2)){
    Table[i,j] <- c(paste(round(T1[i,j],3),"(",round(T1_se[i,j],3),")",sep=""))
  }
}

Table

Table2 <- round(T2,3)
for(i in c(1:5)){
  Table2[i] <- c(paste(round(T2[i],3),"(",round(T2_se[i],3),")",sep=""))
}

Table2

Table3 <- round(T3,3)
for(i in c(1:4)){
  Table3[i] <- c(paste(round(T3[i],3),"(",round(T3_se[i],3),")",sep=""))
}

Table3

require(kableExtra)
kbl(Table3,format = "latex",
    booktabs = T, align = "c")

kbl(cbind(Table[,1],round(T1[1,1]/T1[,1],3),Table[,2],round(T1[1,1]/T1[,2],3)),format = "latex",
    booktabs = T, align = "c")


kbl(cbind(Table2,round(T1[1,1]/PE.pens,3)),format = "latex",
    booktabs = T, align = "c")

kbl(cbind(Table3,round(T1[1,1]/PE.mls,3)),format = "latex",
    booktabs = T, align = "c")

#### plots
require(reshape)
sub1 <- melt(MSE_sub1)
sub2 <- melt(MSE_sub2)
pm   <- melt(MSE_pens)
colnames(MSE_ml)[4] <- c("OLS")
ml   <- melt(MSE_ml)


df1<- data.frame(Model=rep("Shrinkage Methods(via Forward)",dim(sub1)[1]),sub1)
df2<- data.frame(Model=rep("Shrinkage Methods(via Backward)",dim(sub2)[1]),sub2)
df3<- data.frame(Model=rep("Penalized Methods",dim(pm)[1]),pm)
df4<- data.frame(Model=rep("Machine Learnings + OLS",dim(ml)[1]),ml)


df_boxplot<- rbind(df1,df2,df3,df4)



require(viridis)
require(hrbrthemes)
# Plot
plot1 <- ggplot(df_boxplot, aes(x=factor(X2), y=value) ) + 
  #geom_boxplot(alpha = 0.5, show.legend = FALSE, aes(fill=factor(X2))) + 
  geom_boxplot() + 
  #geom_jitter(width=0.25)+
  scale_y_continuous(limits=c(0,15000)) + 
  facet_grid(.~Model,scales = "free") +
  theme_light()+
  theme(strip.text.x = element_text(size=10, color="black", face="bold")) +
  labs(x=" ", y="PE")
plot1

#setwd("/Users/bahadiryuzbasi/Desktop/marwan al-momani/KL_paper_Last_Version")
setEPS()
postscript("box_plot_pe_original_data_1000_rep.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 8, width = 14)
plot1
dev.off()
