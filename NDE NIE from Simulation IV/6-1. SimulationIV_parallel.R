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

####read in simulated datasets by simseed####
for (simseed in 0:199){
  
data=read.csv(paste("sim2_alldata",as.character(simseed+1),".csv",sep=""),header=TRUE)
X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2


####read in estimated parameter  by simseed####


est=read.table(paste("estIVlog",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1)
est = as.data.frame(t(est))


###define parameters
tseq=seq(1,8,by=1)
###number of sampling
B=10000
###number of bootstrap 
M=100
n=length(X)

###get point estimate (save in a vector form)
myest=mymed(X,n,est,tseq,mod=1,B)

###save results
write.csv(myest,paste("estIV_lres",as.character(simseed+1),".csv",sep=""),row.names=FALSE)
}
