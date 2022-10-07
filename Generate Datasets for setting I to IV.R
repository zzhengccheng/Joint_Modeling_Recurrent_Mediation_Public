##Load the package

library(MASS)
library(survival)
library(dplyr)
library(reda)
library(clusterGeneration)
###load in all self defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("myprod.R")
source("datagen_X.R")


################################################################################################
####################################Setting###################################################
n=200
k=1
tau=10
t.jump.m=c(2,4,6)
t.jump.y=c(2,4,6)
lambda.m=c(0.3,0.5,0.6,0.65)
lambda.y=c(0.033,0.084,0.136,0.8)

betax=rep(0.2,k)
betaz=0.35
etaz=0.35
etax=rep(0.15,k)
etam=0.25
delta1=1
sigmav=0.7

#################################Generate 200 datasets from Setting I#########################
delta2=0
for (iter in 1:200){
  set.seed(8888+iter)
  vi=rnorm(n)*sigmav
  data=datagen_X(n,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y)
  write.csv(data$alldata,paste("sim1_alldata",as.character(iter),".csv",sep=""),row.names=F)
}


#################################Generate 200 datasets from Setting II#########################
delta2=-0.5
for (iter in 1:200){
  set.seed(8888+iter)
  vi=rnorm(n)*sigmav
  data=datagen_X(n,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y)
  write.csv(data$alldata,paste("sim2_alldata",as.character(iter),".csv",sep=""),row.names=F)
}


#################################Generate 200 datasets from Setting III#########################
for (iter in 1:200){
  set.seed(8888+iter)
  v_gamma=rgamma(n,shape = 1.5, scale = 1)
  vi=log(v_gamma)
  data=datagen_X(n,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y)
  write.csv(data$alldata,paste("sim3_alldata",as.character(iter),".csv",sep=""),row.names=F)
}

