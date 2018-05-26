### This is the code used to analyze the Cancer dataset in the manuscript


#source("http://www.bioconductor.org/biocLite.R")
#biocLite("GEOquery")

#install.packages(ggplot2)
#install.packages(glmnet)
#install.packages(grplasso)
#install.packages(GSA)

## This function is used to mean impute the missing data

mean.impute <- function(X){
  means <- apply(X,2,function(x){mean(x[which(!is.na(x))])})
  for(i in 1:ncol(X)){
    ind <- which(is.na(X[,i]))
    X[ind,i] <- means[i]
  }
  return(X)
}

## This function checks what proportion of the data is missing

prop.missing <- function(X){
  apply(X,2,function(x){mean(is.na(x))})
}

library(ggplot2)
library(Biobase)
library(GEOquery)
library(glmnet)
library(SGL)
library(grplasso)
library("GSA")

## This grabs the data from bioconductor

gds807 <- getGEO('GDS807', destdir = "~/readingSoft/")

## This preprocesses it

eset <- GDS2eSet(gds807, do.log2 = TRUE)

## This constructs our design matrix and reponse

y <- (as.numeric(pData(eset)$disease.state) < 2)
X <- t(exprs(eset))

## Here we remove all genes with > 50% missingness

prop.m <- prop.missing(X)
remove.ind <- which(prop.m > 0.5)
imp.X <- mean.impute(X[,-remove.ind])
X <- imp.X

## This grabs the gene identifiers

Gene.Identifiers <- Table(gds807)[-remove.ind,2]

## The following code creates the group index using the C1 genesets

filename="C1.gmt"
junk1=GSA.read.gmt(filename)

index <- rep(0,length(Gene.Identifiers))
for(i in 1:277){
  indi <- match(junk1$genesets[[i]],Gene.Identifiers)
  index[indi] <- i
}

Gene.set.info <- junk1  
dim(X)
length(y)
ind.include <- which(index != 0)
genenames <- Gene.Identifiers[ind.include]
X <- X[,ind.include]
membership.index <- rep(0,ncol(X))
for(i in 1:277){ 
  for(j in 1:length(Gene.set.info$genesets[[i]])){
    change.ind <- match(Gene.set.info$genesets[[i]][j],genenames)
    if(!is.na(change.ind)){
      if(membership.index[change.ind] == 0){
        membership.index[change.ind] <- i
      }
    }
  }
}

## In the following section we train our models and test on heldout data

set.seed(100) 

## We choose a training and test set

training.ind <- sample(1:nrow(X), 30)

train.data <- list(x = X[training.ind,], y = y[training.ind])
test.data <- list(x = X[-training.ind,], y = y[-training.ind])

## We standardize the variables for the group lasso

x.gl <- t(t(train.data$x) - apply(train.data$x,2,mean))
x.gl <- t(t(x.gl) / apply(x.gl,2,sd))
x.gl <- cbind(1, x.gl)

## This runs the group lasso code

index.gl <- c(NA, membership.index)
lambda.max.group <- lambdamax(x.gl, as.numeric(train.data$y), index.gl, standardize = FALSE)
lambdas.gl <- exp(seq(from = log(lambda.max.group), to = log(lambda.max.group*0.1), length.out = 100))
fit.gl <- grplasso(x.gl, as.numeric(train.data$y), index.gl, lambda = lambdas.gl, standardize = FALSE)

## We classify held out observations

t.x.gl <- t(t(test.data$x) - apply(test.data$x,2,mean))
t.x.gl <- t(t(t.x.gl) / apply(t.x.gl,2,sd))
t.x.gl <- cbind(1, t.x.gl)

test.pred.GL <- predict(fit.gl, t.x.gl)
test.pred.GL <- exp(test.pred.GL) / (1 + exp(test.pred.GL))

fit <- glmnet(train.data$x, train.data$y, family = "binomial", lambda.min.ratio = 0.1) # 0.01?

test.pred <- predict(fit, test.data$x, type = "class")-1

fitSGL <- SGL(train.data, membership.index, type = "logit", verbose = TRUE, nlam = 100, min.frac = 0.1, alpha = 0.05) 

test.pred.SGL <- matrix(NA, ncol = length(fitSGL$lambdas), nrow = length(test.data$y))

for(i in 1:length(fitSGL$lambdas)){
  test.pred.SGL[,i] <- predict(fitSGL,test.data$x,i)
}

## We see how well each model performed

correct.class <- (test.pred > 0.5) * test.data$y + (test.pred < 0.5)* (1-test.data$y)

correct.class.SGL <- (test.pred.SGL > 0.5) * test.data$y + (test.pred.SGL < 0.5)* (1-test.data$y)

correct.class.GL <- (test.pred.GL > 0.5) * test.data$y + (test.pred.GL < 0.5)* (1-test.data$y)

c.l <- apply(correct.class,2,mean)
c.gl <- apply(correct.class.GL,2,mean)
c.sgl <- apply(correct.class.SGL,2,mean)

max(c.sgl)
max(c.gl)
max(c.l)

best.sgl <- which.max(c.sgl)
best.gl <- which.max(c.gl)
best.l <- which.max(c.l)

best.ind.SGL <- which(fitSGL$beta[,best.sgl] != 0)
best.ind.GL <- which(fitSGL$beta[,best.gl] != 0)
best.ind.l <- which(fit$beta[,best.l] != 0)

Gene.set.info$geneset.names[unique(membership.index[best.ind.SGL])]
Gene.set.info$geneset.names[unique(membership.index[best.ind.GL])]
unique(membership.index[best.ind.l])

c.class <- c(c.gl,c.l,c.sgl)
Method <- c(rep("GL",100),rep("Lasso", 100),rep("SGL",100))
c.lambda <- c(1:100, 1:100,1:100)

## We plot the results

#pdf("cancer.pdf")
dd <- data.frame(Method = Method, x = c.lambda, y = c.class)

ggplot(data = dd, aes(x = x, y = y, group = Method, shape = Method)) + geom_line(aes(linetype=Method), size = 1.5) + scale_y_continuous("Correct Classification Rate") + scale_x_continuous("Lambda Index") + opts(title = "Correct Classification Rate for Cancer Data",legend.text = theme_text(size = 20), plot.title = theme_text(size = 22), axis.title.x = theme_text(size = 20), axis.title.y = theme_text(size = 20, angle = 90)) + scale_linetype_manual(values=c("twodash","dotted","solid")) + opts(legend.key.size = unit(2,"cm"))

#dev.off()
