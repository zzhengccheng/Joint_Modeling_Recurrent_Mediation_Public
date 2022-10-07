
library(MASS)
library(survival)
library(dplyr)
library(reda)
library(clusterGeneration)
###load in self-defined functions#####
source("Lambdainv.R")
source("mysimRec.R")
source("myprod.R")
source("datagen_M.R")
source("datagen.R")
source("S00.R")
source("S01.R")
source("S10.R")
source("S11.R")
source("mymed2.R")


#######################compute true effects using large data#####################################
#################################################################################################
##### TRUE parameter ######
k=1
tau=10
t.jump.m=c(2,4,6)
t.jump.y=c(2,4,6)
n=10000

lambda.m=c(0.3,0.5,0.6,0.65)
lambda.y=c(0.033,0.084,0.136,0.8)
betax=rep(0.2,k)
betaz=0.35
etaz=0.35
etax=rep(0.15,k)
etam=0.25
delta1=1
delta2=-0.5
sigmav=0.7

#####Generate a large dataset#######


set.seed(8888)
vi=rnorm(n)*sigmav
# Generate the corrlation matrix

mu<- rep(0,k)
corMat <- rcorrmatrix(k)

# Generate data X from normal distribution
X <- mvrnorm(n, mu=mu, Sigma=corMat, empirical = TRUE)
data=datagen(n,X,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y)$alldata

X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2
tseq=seq(1,8,by=1)

#####Get the true value####
truevalue=mymed2(n,X,sigmav,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,tseq,B=1000)
truevalue1=c(truevalue$NDE,truevalue$NIE)

#######handle each simulated datasets##########

allest=NULL
allsd=NULL
myest=NULL

for (i in 1:200){
  myBest=NULL
  myest1=read.csv(paste("estIII_lres",as.character(i),".csv",sep=""))
  myest=c(myest1$NDE[-c(9,18)],myest1$NIE[-c(9,18)])
  for (batch in 1:5){
    
    skip_to_next <- FALSE
    tryCatch({myBest_temp=read.csv(paste("BootIII_lres",as.character(i),"batch=",as.character(batch),".csv",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    if (!is.na(sum(myBest_temp))){
      myBest=rbind(myBest,myBest_temp)
    }
  }
  
  
  if (!is.null(myBest)){
    myBest[abs(myBest)>1]<-NA
    myest.sd=apply(myBest,2,sd,na.rm=TRUE)
    allest=rbind(allest,myest)
    allsd=rbind(allsd,myest.sd)
  }
}

bias=apply(allest,2,mean,na.rm=TRUE)-truevalue1
ese=apply(allsd*is.finite(allsd),2,mean,na.rm=TRUE)
meSE=apply(allsd,2,median, na.rm=TRUE)
sd=apply(allest,2,sd,na.rm=TRUE)
cr=apply((abs(allest-rep(1,nrow(allest))%o%truevalue1)/allsd)<1.96,2,mean,na.rm=TRUE)

###Save the result###
write.csv(rbind(truevalue1,bias,ese,sd,cr),"sim3res_log.csv")



