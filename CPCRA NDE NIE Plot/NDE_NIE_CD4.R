
library(MASS)
library(survival)
library(dplyr)
###load in self-defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("datagenQ13_M.R")
source("myprod.R")
source("S.R")
source("mymed_Q13.R")

####check real data###########
####read in data##############

covdata=read.csv("covdata.csv")
covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM
covdata$cd4bl=covdata$CD4BL/100
covdata$prevoi=covdata$PREVOI
############Z is CD4BL##########
Z=covdata[,"cd4bl"]
Q1=unname(quantile(Z,na.rm=TRUE, probs=0.25))
Q3=unname(quantile(Z,na.rm=TRUE, probs=0.75))

tseq=seq(1,15,by=1)

X=as.matrix(covdata[,c("prevoi","gender","stratum","hemobl","trt")])
tau=15
n=nrow(X)

est=read.table(paste("estI_real",".txt",sep=""),header=FALSE,row.names=1)
est = as.data.frame(t(est))

betaz=est$betax5
betax=c(est$betax1,est$betax2,est$betax3,est$betax4,est$betaz)
lambda.m=exp(c(est$log_r1,est$log_r2,est$log_r3,est$log_r4,
               est$log_r5,est$log_r6,est$log_r7,est$log_r8,
               est$log_r9,est$log_r10))
t.jump.m=c(1.067,2.267,3.567,4.767,6.333,7.8,9.067,11.1,13.567)
lambda.y=exp(c(est$log_h1,est$log_h2,est$log_h3,est$log_h4,
               est$log_h5,est$log_h6,est$log_h7,est$log_h8,
               est$log_h9,est$log_h10))
t.jump.y=c(2.4,3.633,5.433,7.3,8.6,10.233,11.133,12.267,14.867)
etaz=est$etax5
etax=c(est$etax1,est$etax2,est$etax3,est$etax4,est$etaz)
etam=est$etam
delta1=est$delta1
delta2=0
sigmav=sqrt(exp(est$log_varc))




B=100000/20
###Load the estimates from real data analysis####
est1=c(est$log_r1,est$log_r2,est$log_r3,est$log_r4,est$log_r5,
       est$log_r6,est$log_r7,est$log_r8,est$log_r9,est$log_r10,
       est$log_h1,est$log_h2,est$log_h3,est$log_h4,est$log_h5,
       est$log_h6,est$log_h7,est$log_h8,est$log_h9,est$log_h10,
       est$betax1,est$betax2,est$betaz,est$betax3,est$betax4,
       est$betax5,est$etax1,est$etax2,est$etaz,est$etax3,est$etax4,
       est$etax5,est$delta1,est$etam,est$log_varc)


##Mediation Analysis##
for (batch in 1:20){
  set.seed(1111)
  Bres1_new=NULL
  Best=est1
  Bsigmav=sqrt(exp(Best[35]))
  Blambda.m=exp(Best[1:10])
  Blambda.y=exp(Best[11:20])
  Bbetaz=Best[26]
  Bbetax=Best[c(21:22,24:25,23)]
  Betaz=Best[32]
  Betax=Best[c(27:28,30:31,29)]
  Bdelta1=Best[33]
  Bdelta2=0
  Betam=Best[34]
  Bres1_new=mymed_Q13(X,n,Z,Q1,Q3,Bsigmav,tau,Bbetaz,Bbetax,Blambda.m,t.jump.m,Betaz,Betax,Betam,Bdelta1,Bdelta2,Blambda.y,t.jump.y,tseq,B,seed=1111+batch*B)
 

 write.csv(Bres1_new,paste("Res_cd4",as.character(1111),"batch=",as.character(batch),".csv",sep=""),row.names=FALSE)
}