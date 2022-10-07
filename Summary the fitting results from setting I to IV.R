
####################################################################
##########################True parameters###########################
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

################Computing log() for lambda.m, lambda.y and sigmav###############

log_lam_m1=log(lambda.m[1])
log_lam_m2=log(lambda.m[2])
log_lam_m3=log(lambda.m[3])
log_lam_m4=log(lambda.m[4])

log_lam_y1=log(lambda.y[1])
log_lam_y2=log(lambda.y[2])
log_lam_y3=log(lambda.y[3])
log_lam_y4=log(lambda.y[4])

log_varc=log(sigmav^2)

##############################################################################
###################################Setting I###################################
True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
            log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
            betax,betaz,etaz,etax,etam,delta1,log_varc)

#####Import the estimation and covariance######


est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  est.next=read.table(paste("estIlog",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  cov=read.table(paste("covIlog",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))
  var.next=t(diag(as.matrix(cov)))
  est=rbind(est,est.next[2,])
  var=rbind(var,var.next)
}
names(est)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","log_varc") 
names(var)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","log_varc") 


est=as.data.frame(sapply(est, as.numeric))

################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("sim1_estbialog",".csv",sep=""),row.names = TRUE)



##############################################################################
###################################Setting II###################################
delta2=-0.5
True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
          log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
          betax,betaz,etaz,etax,etam,delta1,delta2,log_varc)
#####Import the estimation and covariance######

est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  est.next=read.table(paste("estIIlog",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  cov=read.table(paste("covIIlog",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))
  var.next=t(diag(as.matrix(cov)))
  est=rbind(est,est.next[2,])
  var=rbind(var,var.next)
}
names(est)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","delta2","log_varc") 
names(var)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","delta2","log_varc") 


est=as.data.frame(sapply(est, as.numeric))
################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("sim2_estbialog",".csv",sep=""),row.names=TRUE)



##############################################################################
###################################Setting III###################################

True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
          log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
          betax,betaz,etaz,etax,etam,delta1,delta2,log_varc)

#####Import the estimation and covariance######

est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  
  skip_to_next <- FALSE
  tryCatch({cov=read.table(paste("covIIIlog",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))}, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
  var.next=t(diag(as.matrix(cov)))
  var=rbind(var,var.next)
  
  est.next=read.table(paste("estIIIlog",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  est=rbind(est,est.next[2,])
  
}
names(est)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","delta2","log_varc") 
names(var)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","delta2","log_varc") 


est=as.data.frame(sapply(est, as.numeric))

################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)[-16,]
write.csv(my_eva,paste("sim3_estbialog",".csv",sep=""),row.names=TRUE)



##############################################################################
###################################Setting IV#################################

True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
          log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
          betax,betaz,etaz,etax,etam,delta1,log_varc)
#####Import the estimation and covariance######

est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  
  skip_to_next <- FALSE
  tryCatch({cov=read.table(paste("covIVlog",as.character(iter),".txt",sep=""),row.names=1,na.strings=c("."))}, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
  var.next=t(diag(as.matrix(cov)))
  var=rbind(var,var.next)
  
  est.next=read.table(paste("estIVlog",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  est=rbind(est,est.next[2,])
  
}
names(est)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","log_varc") 
names(var)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4","log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4","betax","betaz","etaz","etax","etam","delta1","log_varc") 


est=as.data.frame(sapply(est, as.numeric))
################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE,ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("sim4_estbialog",".csv",sep=""),row.names=TRUE)

