# attempt SVM model generation for reference set integration
# check for installed packages. Install missing as needed
pkgs = c("data.table","e1071","rpart","caret","pROC","randomForest","elasticnet","plyr","kernlab","doParallel","arm")
if(length(new.pkgs <- setdiff(pkgs, rownames(installed.packages())))>0) install.packages(new.pkgs)
rm(pkgs,new.pkgs)
library(pROC)
library(kernlab)
library(gbm)
library(survival)
library(mboost)
library(caret)
library(data.table)
library(plyr)
library(e1071)
library(doParallel)
# depending on computer used, parallelization can be problematic and result
# in R or complete system crashes. Following commented code is in case it is
# desired to try. If so, be sure to stop or unregister the cores, or to
# return processing from parallel to registerDoSEQ()
# unregister <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }

# define by making clusters
# no_cores <- detectCores() - 2
# cl <- makeCluster(4)
## or ##
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# stopCluster(cl)

# define by setting cores and let system handle the details
# registerDoParallel(cores=2)

# stop any clusters/cores
# stopCluster(cl)
# unregister()
# registerDoSEQ()

##############################################################################
# in original study, statements were here included to process the collected set
# of validated interactions and the "Range_floored_master_tab" table of all
# prediction algorithms and corresponding scores. This abbreviated version
# is designed to work with the convenience files already provided in the
# github repository: train_data.txt and test_data.txt

# if(!(exists("master_tab.Fl"))){
#   load("Range_floored_master_tab.RData")
# }
#
# load("all_valid_sets.RData")
# totExptVal <- totExptVal[!is.na(Type)]
#
# # remove the FC variant with fold changes, as classification rather than
# # regression to be performed
# rm(totExptValFC)
#
# # update column types in totExptVal object
# totExptVal[,c("miRID","Type","Source"):=
#              list(as.character(miRID),
#                   factor(Type,levels=c("pos","neg")),
#                   as.factor(Source))]
#
# # create object of master table merged with validated results
# setkey(master_tab.Fl,GeneName,miRID)
# setkey(totExptVal,GeneName,miRID)
# selectML <- totExptVal[master_tab.Fl,nomatch=0]
# # evaluate negative duplicates
# dupes <- selectML[duplicated(selectML)][,c(1:4),with=F]
# temp <- selectML[dupes,nomatch=0][,c(1:4),with=F]
# temp2 <- selectML[!temp]
#
# setkey(temp,GeneName,miRID,Type,Source)
# temp <- unique(temp)
# temp[,"coding":=ifelse(Type=="pos",1,-1)]
# temp[,"check":=sum(coding),by=.(GeneName,miRID)]
#
# # in cases where check==1, there were two pos votes, one neg.
# # remove the sole dissenter. Also, when check==0, take the "VALID" result
# rej <- which(temp[,check]==1&temp[,Type]=="neg")
# rej2 <- which(temp[,check]==0&temp[,Source!="VALID"])
# kee <- setdiff(c(1:nrow(temp)),union(rej,rej2))
# temp <- temp[kee]
# # recalculate check across remaining entries
# temp[,"check":=sum(coding),by=.(GeneName,miRID)]
# confirmed <- temp[abs(check)>1]
# save(confirmed,file="experimentally_confirmed_interactions_multiple_source.RData")
# setkey(temp,GeneName,miRID,Type)
# temp <- unique(temp,fromLast = T)
# temp[,c("coding","check"):=NULL]
# setkey(selectML,GeneName,miRID,Type,Source)
# setkey(temp,GeneName,miRID,Type,Source)
# temp <- selectML[temp,nomatch=0]
# trimmedML <- do.call(rbind,list(temp,temp2))
# rm(temp,temp2,kee,rej,rej2,dupes,confirmed)
# # downsampling to level neg and pos class sizes. Ensure stratified
# # generate table of summaries,
# setkey(trimmedML,Type,Source)
# tabsumm <- trimmedML[,.N,by=key(trimmedML)]
#
# # take median value of negative counts for sampling.
# # will thus take all VALID neg, all OEKD neg (corresponding samplings)
# # and sample 633 from remainder
# midpt <- median(tabsumm[Type=="neg",N])
# set.seed(31415)
# idx <- numeric()
# for (i in levels(trimmedML[,Source])){
#   if(tabsumm[Source==i&Type=="neg",N]<=midpt){
#     count <- tabsumm[Source==i&Type=="neg",N]
#     idx <- union(idx,which(trimmedML[,Source]==i&trimmedML[,Type]=="neg"))
#     idx <- union(idx,sample(which(trimmedML[,Source]==i&trimmedML[,Type]=="pos"),count))
#   } else {
#     count <- midpt
#     idx <- union(idx,sample(which(trimmedML[,Source]==i&trimmedML[,Type]=="pos"),count))
#     idx <- union(idx,sample(which(trimmedML[,Source]==i&trimmedML[,Type]=="neg"),count))
#   }
# }
#
# # for mixed over/under sample scheme, determine target sample size
# # per pos/neg class, per type of source information
# pcCount <- round(nrow(trimmedML)/9.99,0)
# lvSource <- levels(trimmedML$Source)
# DBouts <- list()
# set.seed(123)
# for (i in 1:length(lvSource)){
#   tempmat <- data.frame(trimmedML[trimmedML$Source==lvSource[i]])
#   posidx <- which(tempmat$Type=="pos")
#   negidx <- which(tempmat$Type=="neg")
#   if(length(posidx)>pcCount){
#     newpos <- sample(posidx,pcCount)
#   } else {
#     newpos <- sample(posidx,pcCount,replace=T)
#   }
#   if(length(negidx)>pcCount){
#     newneg <- sample(negidx,pcCount)
#   } else {
#     newneg <- sample(negidx,pcCount,replace=T)
#   }
#   tempmat <- rbind(tempmat[newpos,],tempmat[newneg,])
#   DBouts[[i]] <- tempmat
# }
# trimmedMLmix <- rbind.fill(DBouts)
# set.seed(835)
# idx <- createDataPartition(trimmedMLmix$Type,p=0.7,list=F)
# trainMix <- trimmedMLmix[idx,-c(1,2,4)]
# testMix <- trimmedMLmix[-idx,]
# rm(count,i,midpt,idx,DBouts,tempmat,posidx,negidx,newpos,newneg,lvSource,pcCount)

trainData <- read.table("train_data.txt",header=T)
# models to consider

# gradient boosting machine
#    package: gbm
#    gbm (n.trees, interaction.depth, shrinkage, n.minobsinnode)
#        (50:100)       (1,3,5,7)     (0.01-0.1)      (20)
#
# noosted generalized additive model
#    package: mboost, plyr
#    gamboost ( mstop,      prune
#              (50:500)  ("yes","no")
#
# C5.0
#    package: C50, plyr
#    c50 ( trials,  winnow,       model )
#          (5:100)   (T,F)   ("tree","rules")
#
# Naive Bayes
#    package: klaR
#    nb (fL, usekernel)
#        (1)    (T,F)
#
# Bayesian generalized linear models
#    package: arm
#    bayesglm (none)
#
# Extreme gradient boosting
#    package: xgboost
#    xgbTree (nrounds, max_depth, eta, gamma, colsample_bytree,
#             min_child_weight)
#
# k-Nearest neighbours
#    package: (none)
#    knn (k)
#        (1,3,5,9,15)
#
# Generalized linear model with stepwise feature selection
#    package: MASS
#    glmStepAIC (none)
#
# Support vector machine with linear kernel
#    package: e1071
#    svmLinear2 ( cost )
#
#
# Support vector machine with radial kernel
#    package: kernlab
#    svmRadialCost (C)

set.seed(357)
fitCtrl <- trainControl( ## 10-fold CV
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  summaryFunction = twoClassSummary,
  classProbs = T
)

# NOTE BEFORE EXECUTING
# depending on the system used, the following code may overload system
# resources. Consider saving objects or workspaces at intermediate points
# for example after the training of each model, to preserve work already
# completed in case of system hang up

# several versions of the gradient boosting machine tuning parameter sets
# were tried. Grid4 represented the last set of improvements. Grid5 showed
# poorer performance on resampling (see later evaluation) but all are included
# for completeness
gbmGrid <-  expand.grid(interaction.depth = c(1,2,5,7),
                        n.trees = (1:10)*100,
                        shrinkage = c(0.01,0.1),
                        n.minobsinnode = 20)
gbmGrid2 <-  expand.grid(interaction.depth = c(5,7,9),
                         n.trees = (1:5)*500,
                         shrinkage = c(0.05,0.1,0.15),
                         n.minobsinnode = 20)
gbmGrid3 <-  expand.grid(interaction.depth = 7,
                         n.trees = (1:5)*1000,
                         shrinkage = c(0.1,0.15,0.2),
                         n.minobsinnode = 20)
gbmGrid4 <-  expand.grid(interaction.depth = 7,
                         n.trees = (6:10)*1000,
                         shrinkage = c(0.2,0.3),
                         n.minobsinnode = 20)
gbmGrid5 <-  expand.grid(interaction.depth = c(9,11,13),
                         n.trees = (1:10)*1000,
                         shrinkage = c(0.2),
                         n.minobsinnode = 20)

set.seed(357)
gbmFit <- train(Type ~ ., data = trainData,
                method = "gbm",
                trControl = fitCtrl,
                verbose = FALSE,
                tuneGrid = gbmGrid,
                preProcess = c("center","scale"),
                metric="ROC")
set.seed(357)
gbmFit2 <- train(Type ~ .,data=trainData,
                 method = "gbm",
                 trControl = fitCtrl,
                 verbose = F,
                 tuneGrid = gbmGrid2,
                 preProcess = c("center","scale"),
                 metric="ROC")
set.seed(357)
gbmFit3 <- train(Type ~ .,data=trainData,
                 method = "gbm",
                 trControl = fitCtrl,
                 verbose = F,
                 tuneGrid = gbmGrid3,
                 preProcess = c("center","scale"),
                 metric="ROC")
set.seed(357)
gbmFit4 <- train(Type ~ .,data=trainData,
                 method = "gbm",
                 trControl = fitCtrl,
                 verbose = F,
                 tuneGrid = gbmGrid4,
                 preProcess = c("center","scale"),
                 metric="ROC")
set.seed(357)
gbmFit5 <- train(Type ~ .,data=trainData,
                 method = "gbm",
                 trControl = fitCtrl,
                 verbose = F,
                 tuneGrid = gbmGrid5,
                 preProcess = c("center","scale"),
                 metric="ROC")

# next model:
set.seed(357)
gamFit <- train(Type ~ ., data = trainData,
                method = "gamboost",
                trControl = fitCtrl,
                tuneGrid = expand.grid(mstop=c(50,100,150,500),prune=c("yes","no")),
                preProcess = c("center","scale"),
                metric="ROC")

# next model: Naive Bayes
set.seed(357)
nbFit <- train(Type ~ ., data = trainData,
               method = "nb",
               trControl = fitCtrl,
               tuneGrid = expand.grid(fL=1,usekernel=c(TRUE,FALSE),adjust=c(0.1,1)),
               metric="ROC")
# next model: Bayesian generalized linear models
set.seed(357)
bglmFit <- train(Type ~ .,data = trainData,
                 method="bayesglm",
                 trControl = fitCtrl,
                 metric="ROC",
                 preProcess = c("center","scale"))

# next model: extreme gradient boosting machine
set.seed(357)
xgBFit <- train(Type ~ ., data = trainData,
                 method = "xgbTree",
                 trControl = fitCtrl,
                 tuneGrid = expand.grid(nrounds=c(1:10),
                                        eta=c(0.3,0.7),
                                        min_child_weight=c(5,100),
                                        colsample_bytree=c(0.2,0.4,0.8),
                                        max_depth=c(2,5),
                                        gamma=0.7),
                 metric="ROC",
                 preProcess = c("center","scale"))

# next model: k-nearest neighbours
set.seed(357)
knnFit <- train(Type ~ ., data = trainData,
                method = "knn",
                trControl = fitCtrl,
                tuneGrid = expand.grid(k=c(1,3,5,10,15,20)),
                metric="ROC")

# next model: Generalized linear model with stepwise feature selection
set.seed(357)
glmSAICFit <- train(Type ~ .,data = trainData,
                    method="glmStepAIC",
                    tuneLength=5,
                    trControl = fitCtrl,
                    metric="ROC",
                    preProcess = c("center","scale"),
                    verbose=F)

# next model: C5.0
set.seed(357)
C50Fit <- train(Type ~ .,data=trainData,
                method="C5.0",
                trControl = fitCtrl,
                tuneGrid = expand.grid(trials=c(5,10,50,100),winnow=c(T,F),model=c("tree","rules")),
                metric="ROC")

# next model: Support vector machine with linear kernel
set.seed(357)
svmL2Fit <- train(Type ~ ., data = trainData,
                  method = "svmLinear2",
                  trControl = fitCtrl,
                  tuneGrid = expand.grid(cost=c(1,10,100)),
                  preProcess = c("center","scale"),
                  metric="ROC")


# next model: Support vector machine with radial kernel, cost-sensitive
set.seed(357)
svmRCFit <- train(Type ~ ., data = trainData,
                  method = "svmRadialCost",
                  trControl = fitCtrl,
                  tuneLength = 3,
                  preProcess = c("center","scale"),
                  metric="ROC")

# use caret in-built resampling function to compare model fitting results
# (resampling of the fitting process. This comparison is only with data
# seen during the training of each model)
models <- list(GBM = gbmFit4,
               stepAIC = glmSAICFit,
               BayesGLM = bglmFit,
               NaiveBayes = nbFit,
               knn = knnFit,
               GAM = gamFit,
               C5.0 = C50Fit,
               xGBM = xgBFit,
               svmLinear = svmL2Fit,
               svmRadial = svmRCFit)
resampMix <- resamples(models)
# generate lattice plot of resampling results
bwplot(resampMix,layout=c(3,1))

# load test data for assessment on unseen data (accuracy and kappa)
testMix <- read.table("test_data_inputs.txt",header=T)
testresults <- lapply(models,
                      function(x){
                        confusionMatrix(predict.train(x,testMix[,-c(1:4)]),
                                        testMix$Type)})
testresults.mat <- t(sapply(testresults,function(x){x$overall[c("Accuracy","Kappa")]},USE.NAMES = F))
print(testresults.mat)

# best models gbm and knn. Use to predict all interactions
# from master_tab.Fl
gbmFit <- models[["GBM"]]
knnFit <- models[["knn"]]
rm(models,testresults,testMix)
gc()
load("Range_floored_master_tab.RData")
# generate class probabilities for GBM and knn models
predgbm <- predict(gbmFit,master_tab.Fl[,c(3:9),with=F],type="prob")
predknn <- predict(knnFit,master_tab.Fl[,c(3:9),with=F],type="prob")

mGgbm <- do.call(cbind,list(master_tab.Fl[,.(GeneName,miRID)],predgbm))
mGknn <- do.call(cbind,list(master_tab.Fl[,.(GeneName,miRID)],predknn))
# clear the large master table object
rm(master_tab.Fl)
# merge into single dataset
# the gbmpos/knnpos columns indicate probability of interaction's
# belonging to the positive class, gbmneg/knnneg of negative class memebership
names(mGgbm)[3:4] <- c("gbmpos","gbmneg")
names(mGknn)[3:4] <- c("knnpos","knnneg")
setkey(mGgbm,GeneName,miRID)
setkey(mGknn,GeneName,miRID)
mGamalg <- merge(mGgbm,mGknn)
# clear intermediate objects
rm(mGgbm,mGknn)
gc()

# effectively take average of the class probabilities from knn and GBM
# models. Shift values such that scores will range from -1 to +1
gkM <- mGamalg[,.(GeneName,miRID,Score=(gbmpos+knnpos)-1)]
