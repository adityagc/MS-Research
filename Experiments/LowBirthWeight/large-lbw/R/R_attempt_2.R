library(grpreg)
library(caret)
library(ggplot2)
library(caTools)
# install.packages("fscaret")
library(fscaret)
data(Birthwt)
hist(Birthwt$bwt)
X <- Birthwt$X
Y <- Birthwt$bwt
set.seed(0)
partition <- createDataPartition(Y, times = 1, p = 0.7, list = FALSE)
Ytrain <- Y[partition]
Ytest <- Y[-partition]
Xtrain <- X[partition, 1:16]
Xtest <- X[-partition, 1:16]
groups <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
grpmodel <- grpreg(Xtrain, Ytrain, group=groups, penalty =  "grLasso", returnY = TRUE,lambda =  0.03)
Ypred <- predict(grpmodel, Xtest)
# metrics = accuracy(Ytest, Ypred, threshold = 0.5)
confusionMatrix(Ypred, Ytest)
MSE(Ypred, Ytest)
# plot(Y_fit, log.l = FALSE)
# summary(Y_fit)
# weights <- Y_fit$fit$beta
# # xtab <- table(Y_fit$Y, Y)
# length(Y)
# length(Y_fit$Y)
# # weights
# write.csv(weights, "weights.csv")