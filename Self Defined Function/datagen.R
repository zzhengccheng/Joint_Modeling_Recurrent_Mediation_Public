####data gen with fixed X
datagen<-function(n,X,vi,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,cen=TRUE){
  
  
  # Generate binary trt
  Z <- rbinom(n,1,0.5)
  p=k
  
  ####Simulate M^z(t)
  ID=1:n
  beta.all=c(betaz,betax,1)
  risk0=c(exp(cbind(0,X,vi)%*%beta.all))
  risk1=c(exp(cbind(1,X,vi)%*%beta.all))
  ####generate survival time until all individual has follow-up at lest time tau
  M0=mysimRec(t.jump.m,lambda.m,tau,risk0,recurrent=TRUE)
  M1=mysimRec(t.jump.m,lambda.m,tau,risk1,recurrent=TRUE)
  ####Simulate T^{zM^z'} one at a time
  etam=c(etam,delta2)
  T00=T01=T10=T11=rep(0,n)
  risk0.y=c(exp(cbind(0,X,vi)%*%c(etaz,etax,delta1)))
  risk1.y=c(exp(cbind(1,X,vi)%*%c(etaz,etax,delta1)))
  for (i in 1:n){
    tmpM0=M0[which(M0$ID==i),]
    tmpM1=M1[which(M1$ID==i),]
    lambda.risk00=exp((0:length(tmpM0$stoptime))*(etam[1]+etam[2]*vi[i]))
    t.jump.risk00=tmpM0$stoptime
    lambda.risk10=exp((0:length(tmpM0$stoptime))*(etam[1]+etam[2]*vi[i]))
    t.jump.risk10=tmpM0$stoptime	
    lambda.risk01=exp((0:length(tmpM1$stoptime))*(etam[1]+etam[2]*vi[i]))
    t.jump.risk01=tmpM1$stoptime
    lambda.risk11=exp((0:length(tmpM1$stoptime))*(etam[1]+etam[2]*vi[i]))
    t.jump.risk11=tmpM1$stoptime	
    my00=myprod(t.jump.y,lambda.y,t.jump.risk00,lambda.risk00)
    my10=myprod(t.jump.y,lambda.y,t.jump.risk10,lambda.risk10)	
    my01=myprod(t.jump.y,lambda.y,t.jump.risk01,lambda.risk01)
    my11=myprod(t.jump.y,lambda.y,t.jump.risk11,lambda.risk11)	
    T00[i]=mysimRec(my00[[1]],my00[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T10[i]=mysimRec(my10[[1]],my10[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
    T01[i]=mysimRec(my01[[1]],my01[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T11[i]=mysimRec(my11[[1]],my11[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
  }
  ####Simulate Z
  ID1=ID[which(Z==1)]
  ID0=ID[which(Z==0)]
  Mdat=rbind(M0[which(M0$ID%in%ID0),],M1[which(M1$ID%in%ID1),])
  Mdat=Mdat[which(Mdat$stoptime<=tau),]
  Mdat=Mdat[order(Mdat$ID,Mdat$stoptime),]
  T=T00*(1-Z)+T11*Z
  C=tau-0.01
  if (cen){
    C=runif(n,tau/5,tau)
  }
  Y=pmin(C,T)
  D=as.numeric(T<=C)
  ydata=data.frame(ID,Y,D,Z,X)
  names(ydata)=c("ID","Y","D",paste("X",as.character(1:(p+1)),sep=""))
  mdata=Mdat
  adata=merge(ydata,mdata,by="ID")
  adata=adata[which(adata$stoptime<adata$Y),]
  adata$D=adata$event
  adata$Y=adata$stoptime
  adata=adata[,c("ID","Y","D",paste("X",as.character(1:(p+1)),sep=""))]
  ydata$D=ydata$D*2
  alldata=rbind(adata,ydata)
  alldata=alldata[order(alldata$ID,alldata$Y),]
  ###first is Z, all others are X
  names(alldata)=c("ID","stoptime","event",paste("X",as.character(1:(p+1)),sep=""))
  ydata$vi=vi
  data=list(alldata=alldata,edat=data.frame(ID,T00,T01,T10,T11),mdat=merge(data.frame(Mdat),ydata[,-c(2,3)],by="ID"))
}
