
library(MASS)
library(survival)
library(dplyr)
###load in self-defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("datagen_M.R")
source("S00.R")
source("S01.R")
source("S10.R")
source("S11.R")
source("mymed.R")

###setting running parameters###


for (simseed in 0:999){

bootbatch=simseed%%5+1
simseed=simseed%/%5

####read in simulated datasets by simseed####

data=read.csv(paste("sim1_alldata",as.character(simseed+1),".csv",sep=""),header=TRUE)
X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2
n=length(X)

####read in estimated  parameters by simseed####

est=read.table(paste("estIlog",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1)
est = as.data.frame(t(est))
est_m=c(est$log_lam_m1,est$log_lam_m2,est$log_lam_m3,est$log_lam_m4,
        est$log_lam_y1,est$log_lam_y2,est$log_lam_y3,
        est$log_lam_y4,est$betax,est$betaz,est$etaz,
        est$etax,est$etam,est$delta1,est$log_varc)

####read in estimated variance by simseed####

estvar=read.table(paste("covIlog",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1,na.strings = ".")

###define parameters
tseq=seq(1,8,by=1)
###number of sampling 
B=10000
###number of bootstrap 
M=20

###get point estimate (save in a vector form)

myest=mymed(X,n,est,tseq,mod=1,B)

myBest=matrix(data=NA,nrow=M,ncol=length(myest$NDE)*2)


###compute variance via Bootstrap
if (!is.na(sum(estvar))){
for (m in 1:M){
  
	set.seed(8888+(bootbatch-1)*M+m)
	estb=mvrnorm(1,mu=est_m,Sigma=estvar)
	est=as.data.frame(t(estb))
	mymed_m=mymed(X,n,est,tseq,mod=1,B)
	myBest[m,]=c(mymed_m$NDE,mymed_m$NIE)
	
}
}
###save results

write.csv(myBest,paste("BootI_lres",as.character(simseed+1),"batch=", as.character(bootbatch),".csv",sep=""),row.names=FALSE)
}


