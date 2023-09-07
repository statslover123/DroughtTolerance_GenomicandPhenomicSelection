library('dplyr', lib.loc='/data/lab/pumphrey/RPack/')
library('glmnet')
library('tidyverse')
library('caret', lib.loc='/data/lab/pumphrey/RPack/')
library('tibble')
setwd("/data/lab/pumphrey/RFiles/")
MyY<-read.delim('SpillmanHardandSoftLaergeBlock.txt')
MyG<-read.delim('GAPIT.Genotype.Numerical.txt')
yield <-MyY[c('Taxa','TW')]
yield2 <-MyY[c('Taxa','TW','Yield','NDVI','NDRE1','NWI1','NWI2')]
yield3 <-MyY[c('Taxa','TW','NDVI','NDRE1','NWI1','NWI2')]

names(yield)[1] <- 'taxa'
names(yield2)[1] <- 'taxa'
names(yield3)[1] <- 'taxa'
amp_up_models <- function(){
  library(parallel)
  library(doParallel)
  no_cores <- 4
  #Leave one core available for Operating system
  cluster <- makePSOCKcluster(no_cores)
  registerDoParallel(cluster)
  cat("Model amped and ready to go with:", no_cores, "cores. \n")
}

amp_up_models()  

table2.df <- dplyr::inner_join(yield, MyG, by= 'taxa')
table3.df = subset(table2.df, select = -c(taxa))

table22.df <- dplyr::inner_join(yield2, MyG, by= 'taxa')
table33.df = subset(table22.df, select = -c(taxa))
table33.df <- table33.df %>%
  select(TW, everything())

table222.df <- yield3
table333.df = subset(table222.df, select = -c(taxa))

preObj <- preProcess(table3.df[, -1], method=c("center", "scale"))
newData <- predict(preObj, table3.df[, -1])
newData$TW <- table3.df$TW
newData<-na.omit(newData)


preObj2 <- preProcess(table33.df[, -1], method=c("center", "scale"))
newData2 <- predict(preObj2, table33.df[, -1])
newData2$TW <- table33.df$TW
newData2<-na.omit(newData2)

preObj3 <- preProcess(table333.df[, -1], method=c("center", "scale"))
newData3 <- predict(preObj3, table333.df[, -1])
newData3$TW <- table333.df$TW
newData3<-na.omit(newData3)

lambda.grid <- seq(0, 500, length = 500)
alpha.grid <- seq(0, 1, length = 10)
srchGrd = expand.grid(.alpha = alpha.grid, .lambda = lambda.grid)

cv_5 = trainControl(method = 'repeatedcv', number = 5, repeats=20,savePredictions = 'final',allowParallel = TRUE)

hit_elnet = train(
  TW ~ ., data = newData,
  method = 'glmnet',
  trControl = cv_5,
   tuneGrid = srchGrd
)

hit_elnet2 = train(
  TW ~ ., data = newData2,
  method = 'glmnet',
  trControl = cv_5,
   tuneGrid = srchGrd
)


hit_elnet3 = train(
  TW ~ ., data = newData3,
  method = 'glmnet',
  trControl = cv_5,
   tuneGrid = srchGrd
)

a<-hit_elnet[["resample"]]
a1<-hit_elnet2[["resample"]]
a2<-hit_elnet3[["resample"]]
a <- a %>%
  add_column(Model = "only marker TW")
a1 <- a1 %>%
  add_column(Model = "both marker and drone TW")
a2 <- a2 %>%
  add_column(Model = "only drone TW")
a3<-bind_rows(a, a1,a2)

write.csv(a3, "/data/lab/pumphrey/RScript/ElasticNetWork/TW5foldCV.csv", row.names=FALSE)

pred<-hit_elnet$pred
pred1<-hit_elnet2$pred
pred2<-hit_elnet3$pred


pred <- pred %>%
  add_column(Model = "only marker TW")
pred1 <- pred1 %>%
  add_column(Model = "both marker and drone TW")
pred2 <- pred2 %>%
  add_column(Model = "only drone TW")
pred3<-bind_rows(pred,pred1,pred2)

write.csv(pred3, "/data/lab/pumphrey/RScript/ElasticNetWork/TW5foldCVpred.csv", row.names=FALSE)

tune<-hit_elnet$bestTune
tune2<-hit_elnet2$bestTune
tune3<-hit_elnet3$bestTune

print('tune marker')
tune



print('tune both')
tune2


print('tune drone')
tune3

#fit <- glmnet(TW ~ ., data = newData)
#fit2 <- glmnet(TW ~ ., data = newData2)
#fit3 <- glmnet(TW ~ ., data = newData3)


#coefs<-as.data.frame(as.matrix(coef(fit, s = log(tune$lambda), alpha=tune$alpha)))
#coefs2<-as.data.frame(as.matrix(coef(fit2, s = log(tune2$lambda), alpha=tune2$alpha)))
#coefs3<-as.data.frame(as.matrix(coef(fit3, s = log(tune3$lambda), alpha=tune3$alpha)))
#coefs <- tibble::rownames_to_column(coefs, "Predictor")
#coefs2 <- tibble::rownames_to_column(coefs, "Predictor")
#coefs3 <- tibble::rownames_to_column(coefs2, "Predictor")



coefs<-as.data.frame(as.matrix(coef(hit_elnet$finalModel, hit_elnet$finalModel$lambdaOpt)))
coefs2<-as.data.frame(as.matrix(coef(hit_elnet2$finalModel, hit_elnet2$finalModel$lambdaOpt)))
coefs3<-as.data.frame(as.matrix(coef(hit_elnet3$finalModel, hit_elnet3$finalModel$lambdaOpt)))
coefs <- tibble::rownames_to_column(coefs, "Predictor")
coefs2 <- tibble::rownames_to_column(coefs2, "Predictor")
coefs3 <- tibble::rownames_to_column(coefs3, "Predictor")

coefs <- coefs %>%
  add_column(Model = "only marker TW")
coefs2 <- coefs2 %>%
  add_column(Model = "both marker and drone TW")
coefs3 <- coefs3 %>%
  add_column(Model = "only drone TW")

write.csv(coefs, "/data/lab/pumphrey/RScript/ElasticNetWork/TWonlymarkercoefs.csv", row.names=FALSE)
write.csv(coefs2, "/data/lab/pumphrey/RScript/ElasticNetWork/TWmarkeranddronecoefs.csv", row.names=FALSE)
write.csv(coefs3, "/data/lab/pumphrey/RScript/ElasticNetWork/TWonlydronecoefs.csv", row.names=FALSE)