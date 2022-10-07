S10<-function(Xi,vi,R1i,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,tseq){
  Cb0=c(lambda.y[1],lambda.y[2:length(lambda.y)]-lambda.y[1:(length(lambda.y)-1)])
  sss=length(tseq)
  S10_seq=rep(0,sss)
  for (seq in 1:sss){
    
    j_prime=nrow(R1i)
    sum_j=matrix(0,j_prime+1,4)
    
    for (j in 0:3){
      
      if (tseq[seq]-c(0,t.jump.y)[j+1]>0){sum_j[1,j+1]=Cb0[j+1]*(tseq[seq]-c(0,t.jump.y)[j+1])}
      if (tseq[seq]-c(0,t.jump.y)[j+1]<=0){sum_j[1,j+1]=0}
      
      if (j_prime!=0){
        for (j_p in 1:j_prime){
          T=tseq[seq]-max(c(0,t.jump.y)[j+1],R1i$stoptime[j_p],na.rm=TRUE)
          Cbj=exp((etam+delta2*vi)*(j_p-1))*(exp(etam+delta2*vi)-1)
          
          if (T <= 0){sum_j[j_p + 1,j + 1]=0}
          if (T > 0){sum_j[j_p + 1,j + 1] = Cb0[j+1] * Cbj * T}
          
        }
      }
      sum_2=sum(sum_j)
      
    }
    
    S=exp(-exp(etax*Xi+delta1*vi)*sum_2)
    
    S10_seq[seq]=S
  }
  
  return(S10_seq)
}