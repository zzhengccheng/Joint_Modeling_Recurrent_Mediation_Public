# Compute the True NDE/NIE
mymed2<-function(n,X,sigmav,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,tseq,B){
  NDEmat=NIEmat=matrix(data=NA,nrow=B,ncol=length(tseq))
  
  for (b in 1:B){
    set.seed(1111+b)
    
    vi=rnorm(n)*sigmav
    edata=datagen(n,X,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y)[[2]]
    
    NDE0=apply(outer(edata$T10,tseq,">"),2,mean)-apply(outer(edata$T00,tseq,">"),2,mean)
    NDE1=apply(outer(edata$T11,tseq,">"),2,mean)-apply(outer(edata$T01,tseq,">"),2,mean)
    NIE0=apply(outer(edata$T01,tseq,">"),2,mean)-apply(outer(edata$T00,tseq,">"),2,mean)
    NIE1=apply(outer(edata$T11,tseq,">"),2,mean)-apply(outer(edata$T10,tseq,">"),2,mean)
    NDEmat[b,]=(NDE0+NDE1)/2
    NIEmat[b,]=(NIE0+NIE1)/2
  }
  TEmat=NDEmat+NIEmat
  list(NDE=apply(NDEmat,2,mean,na.rm=TRUE),NIE=apply(NIEmat,2,mean,na.rm=TRUE),TE=apply(TEmat,2,mean,na.rm=TRUE),NDE.SE=apply(NDEmat,2,sd,na.rm=TRUE)/sqrt(B),NIE.SE=apply(NIEmat,2,sd,na.rm=TRUE)/sqrt(B),TE.SE=apply(TEmat,2,sd,na.rm=TRUE)/sqrt(B))
}


