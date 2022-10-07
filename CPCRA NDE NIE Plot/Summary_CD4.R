library(sas7bdat)
library(survival)

MM=20
tseq=seq(0,15*30,by=30)

#########################################################
##########First summary the real estimates of NIE/NDE####
#########################################################
res_NDE=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_NIE=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_TE=matrix(data=NA,nrow=MM,ncol=length(tseq))
res=NULL
for (mm in 1:MM){
  res=read.csv(paste("Res_cd4",as.character(1111),"batch=",as.character(mm-1),".csv",sep=""))
  res_NDE[mm,]=c(0,res$NDE)
  res_NIE[mm,]=c(0,res$NIE)
  res_TE[mm,]=c(0,res$TE)
  
}

R_NDE=apply(res_NDE,2,mean,na.rm=TRUE)
R_NIE=apply(res_NIE,2,mean,na.rm=TRUE)
R_TE=apply(res_TE,2,mean,na.rm=TRUE)


################################################
##########Summary the bootstrap of NIE/NDE######
################################################

Ba_tmp=20
simseed=1000
tmpres_NDE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpres_NIE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpres_TE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
Ares_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
Ares_NIE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
Ares_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_NIE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))


Bres_tmp=NULL
for (sim in 1:simseed){
  for (mm in 1:Ba_tmp){
    Bres_tmp=read.csv(paste("Bres_cd4",as.character(1111+sim-1),"batch=",as.character(mm),".csv",sep=""))
    tmpres_NDE[mm,]=c(0,Bres_tmp$NDE)
    tmpres_NIE[mm,]=c(0,Bres_tmp$NIE)
    tmpres_TE[mm,]=c(0,Bres_tmp$TE)
    
  }
  Ares_NDE[sim,]=apply(tmpres_NDE,2,mean,na.rm=TRUE)
  Ares_NIE[sim,]=apply(tmpres_NIE,2,mean,na.rm=TRUE)
  Ares_TE[sim,]=apply(tmpres_TE,2,mean,na.rm=TRUE)
  SEres_NDE[sim,]=apply(tmpres_NDE,2,sd,na.rm=TRUE)
  SEres_NIE[sim,]=apply(tmpres_NIE,2,sd,na.rm=TRUE)
  SEres_TE[sim,]=apply(tmpres_TE,2,sd,na.rm=TRUE)
  
}


SD_NDE=apply(Ares_NDE,2,sd,na.rm=TRUE)
SD_NIE=apply(Ares_NIE,2,sd,na.rm=TRUE)
SD_TE=apply(Ares_TE,2,sd,na.rm=TRUE)

SENDE=apply(SEres_NDE,2,mean,na.rm=TRUE)
SENIE=apply(SEres_NIE,2,mean,na.rm=TRUE)
SETE=apply(SEres_TE,2,mean,na.rm=TRUE)

QL_NDE=apply(Ares_NDE,2,quantile, probs=0.025,na.rm=TRUE)
QU_NDE=apply(Ares_NDE,2,quantile, probs=0.975,na.rm=TRUE)

QL_NIE=apply(Ares_NIE,2,quantile, probs=0.025,na.rm=TRUE)
QU_NIE=apply(Ares_NIE,2,quantile, probs=0.975,na.rm=TRUE)

QL_TE=apply(Ares_TE,2,quantile, probs=0.025,na.rm=TRUE)
QU_TE=apply(Ares_TE,2,quantile, probs=0.975,na.rm=TRUE)




######################################################
##################Alternative methods#################

covdata=read.sas7bdat("covar.sas7bdat")
covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM
covdata$prevoi=covdata$PREVOI
Z=covdata[,"CD4BL"]/100
Q1=unname(quantile(Z,na.rm=TRUE, probs=0.25))
Q3=unname(quantile(Z,na.rm=TRUE, probs=0.75))

surv=NULL
surv=read.sas7bdat("surv.sas7bdat")
surv$event=surv$event*(1/2)
surv$cd4bl=Z


###############Cox###############
surv=cbind(surv,covdata)
cox_fit <- coxph(Surv(fup*30,event)~cd4bl+trt+gender+hemobl+stratum+prevoi, data=surv)

data0=data1=surv
data0$cd4bl=Q1
data1$cd4bl=Q3

fit2_1<-survfit(cox_fit,newdata=data1)
fit2_0<-survfit(cox_fit,newdata=data0)
Cox_S1<-summary(fit2_1, time=tseq)
Cox_S0<-summary(fit2_0, time=tseq)
fit2<-apply(Cox_S1[["surv"]]-Cox_S0[["surv"]],1,mean,na.rm=TRUE)



###Plot the NDE/NIE/TE ##
pdf("Effects_mod1_cd4_95quantile_final.pdf")
plot(c(R_NDE,R_NIE,R_TE)~c(tseq,tseq,tseq),type="n",ylim=c(-0.10,0.25),ylab="Effect (Difference in Survival Probability)"
     ,xlab="Time (days)",pch=24,cex.lab=1.2)
lines(R_NDE~tseq,col="red",lty=1,lwd=2)
lines(I(QL_NDE)~tseq,col="red",lty=2,lwd=2)
lines(I(QU_NDE)~tseq,col="red",lty=2,lwd=2)

lines(R_NIE~tseq,col="blue",lty=1,lwd=2)
lines(I(QL_NIE)~tseq,col="blue",lty=2,lwd=2)
lines(I(QU_NIE)~tseq,col="blue",lty=2,lwd=2)

lines(R_TE~tseq,col="black",lty=1,lwd=2)
lines(I(QL_TE)~tseq,col="black",lty=2,lwd=2)
lines(I(QU_TE)~tseq,col="black",lty=2,lwd=2)

lines(fit2~tseq,col="green",lty=1,lwd=2)

legend("bottomleft",c("NDE","NIE","TE","Cox"),col=c("red","blue","black","green"),lty=1,lwd=2)
abline(h=0)
dev.off()



