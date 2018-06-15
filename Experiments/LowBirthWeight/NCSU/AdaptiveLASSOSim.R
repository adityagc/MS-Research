##############################################################
##File used for simulation of adaptive lasso and lasso 
#############################################################

#Simulation is from original LASSO paper

set.seed(1314781)
#Need packages for the functions and data creation
library(parcor)
library(mvtnorm)

#number of data sets
M<-500

#True beta vector
beta<-c(3,1.5,0,0,2,0,0,0)
#Define number of predictors
p<-length(beta)
#Create covariance matrix for the X's
a<-diag(p)
sig<-0.8^abs(row(a)-col(a))
n<-20


#################################################
#Create data and standardize it, sample size is 20 for each data set, true sd=3
X<-list()
Y<-list()

for (i in 1:M){
X[[i]]<-rmvnorm(n=n,mean=c(rep(0,p)),sigma=sig)
#standardize X's
for (j in 1:p){
    X[[i]][,j]<-(X[[i]][,j]-mean(X[[i]][,j]))/sd(X[[i]][,j])
}

Y[[i]]<-X[[i]]%*%beta
Y[[i]]<-Y[[i]]+rnorm(n,sd=3)
Y[[i]]<-Y[[i]]-mean(Y[[i]])
}


########################################################################
#Do analysis on each data set

#Create dummy matrices needed to hold the solutions
OLS<-LASSO<-ADLASSO<-matrix(c(0),ncol=p,nrow=M)

for(i in 1:M){
    fits<-adalasso(X=X[[i]],y=Y[[i]],k=5)
    LASSO[i,]<-fits$coefficients.lasso
    ADLASSO[i,]<-fits$coefficients.adalasso
    OLS[i,]<-lm(Y[[i]]~X[[i]]-1)$coef
    
    print(i)
    flush.console()
}

##################################################################
#Now do analysis of how well each method did in terms of mean/median Prediction error and selecting true model
#Also look at how often method selected true model

#Make dummy vectors to hold Model error values
MEOLS<-MELASSO<-MEADLASSO<-c(rep(0,M))

#Now find Model errors
for (i in 1:M){
MEOLS[i]<-t(OLS[i,]-beta)%*%sig%*%(OLS[i,]-beta)
MELASSO[i]<-t(LASSO[i,]-beta)%*%sig%*%(LASSO[i,]-beta)
MEADLASSO[i]<-t(ADLASSO[i,]-beta)%*%sig%*%(ADLASSO[i,]-beta)
}

#Check on model structure
Oracle<-matrix(c(1,1,0,0,1,0,0,0),ncol=p,nrow=M,byrow=TRUE)
OracleLASSO<-sum(apply(FUN=sum,X=((LASSO>0)==Oracle),MARGIN=1)==8)/M
OracleADLASSO<-sum(apply(FUN=sum,X=((ADLASSO>0)==Oracle),MARGIN=1)==8)/M

#mean ME's/median ME's
matrix(c(round(mean(MEOLS),3),round(mean(MELASSO),3),round(mean(MEADLASSO),3),   
         round(sd(MEOLS)/sqrt(M),3),round(sd(MELASSO)/sqrt(M),3),round(sd(MEADLASSO)/sqrt(M),3),
         round(median(MEOLS),3),round(median(MELASSO),3),round(median(MEADLASSO),3),
         0,OracleLASSO,OracleADLASSO),nrow=3,ncol=4,
         dimnames=list(c("OLS","LASSO","ADLASSO"),c("Mean ME","SE","Median ME","Oracle")),byrow=FALSE)

       


###################################################################################
###################################################################################
#Now increase the sample size, for simplicity, just rerun the top with n=100





