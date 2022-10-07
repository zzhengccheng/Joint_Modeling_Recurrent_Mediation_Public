####Function used to compute the NDE/NIE with Z containing two categories Q1 and Q3
mymed_Q13<-function(X,n,Z,Q1,Q3,sigmav,tau,betaz,betax,lambda.m,t.jump.m,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,tseq,B,seed){
  
  #Generate the M_Q1/M_Q3
  avgres=matrix(data=0,nrow=n,ncol=4*length(tseq))  
  avgressq= matrix(data=0,nrow=n,ncol=4*length(tseq))
  
  for (b in 1:B){
    
    set.seed(seed+b)
    vi=rnorm(n,0,sigmav)
    M_T=datagenQ13_M(X,Q1,Q3,vi,betaz,betax,tau,lambda.m,t.jump.m)
    R_Q1=M_T[which((M_T$M==Q1)&(M_T$event==1)),][,-c(3,4)]
    R_Q3=M_T[which((M_T$M==Q3)&(M_T$event==1)),][,-c(3,4)]
    tmpres=avgres
    
    
    for (i in 1:n){
      Xi=X[i,]
      R_Q1i=R_Q1[which(R_Q1$ID==i),]
      R_Q3i=R_Q3[which(R_Q3$ID==i),]
      
      
      S00b=S(Xi,Q1,vi[i],R_Q1i,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,tseq)#modified
      S01b=S(Xi,Q3,vi[i],R_Q1i,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,tseq)#modified
      S10b=S(Xi,Q1,vi[i],R_Q3i,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,tseq)#modified
      S11b=S(Xi,Q3,vi[i],R_Q3i,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,tseq)#modified
      tmpres[i,]= c(S00b,S01b,S10b,S11b)
      
    }
    
    tmpressq=tmpres^2
    
    avgres=avgres+(tmpres-avgres)/b
    avgressq=avgressq+(tmpressq-avgressq)/b
    
  }
  avgsd=sqrt(avgressq-avgres^2)
  
  avg_S=apply(avgres,2,mean,na.rm=TRUE)
  sd_S=apply(avgres,2,sd,na.rm=TRUE)
  
  S00_n=avg_S[1:length(tseq)]
  S01_n=avg_S[(length(tseq)+1):(length(tseq)*2)]
  S10_n=avg_S[(length(tseq)*2+1):(length(tseq)*3)]
  S11_n=avg_S[(length(tseq)*3+1):(length(tseq)*4)]
  
  S00_sd=sd_S[1:length(tseq)]
  S01_sd=sd_S[(length(tseq)+1):(length(tseq)*2)]
  S10_sd=sd_S[(length(tseq)*2+1):(length(tseq)*3)]
  S11_sd=sd_S[(length(tseq)*3+1):(length(tseq)*4)]
  
  
  NDE01=S01_n-S00_n
  NDE10=S10_n-S11_n
  NIE01=S11_n-S01_n
  NIE10=S00_n-S10_n
  
  NDE=(NDE01-NDE10)/2
  NIE=(NIE01-NIE10)/2
  TE=NDE+NIE
  
  return(as.data.frame(list(NDE=NDE,NIE=NIE,TE=TE,S00=S00_n,S01=S01_n,S10=S10_n,S11=S11_n,S00_sd=S00_sd,S01_sd=S01_sd,S10_sd=S10_sd,S11_sd=S11_sd)))
  
}