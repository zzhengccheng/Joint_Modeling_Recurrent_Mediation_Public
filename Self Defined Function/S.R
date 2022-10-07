###Compute the survival functions
S<-function(Xi,Q,vi,Ri,etaz,etax,etam,delta1,delta2,lambda.y,t.jump.y,tseq){
  Cb0=c(lambda.y[1],lambda.y[2:length(lambda.y)]-lambda.y[1:(length(lambda.y)-1)])
  sss=length(tseq)
  S_seq=rep(0,sss)
  J=length(lambda.y)
  for (seq in 1:sss){
    
    j_prime=nrow(Ri)
    sum_j=matrix(0,j_prime+1,J)
    
    for (j in 0:(J-1)){
      
      if (tseq[seq]-c(0,t.jump.y)[j+1]>0){sum_j[1,j+1]=Cb0[j+1]*(tseq[seq]-c(0,t.jump.y)[j+1])}
      if (tseq[seq]-c(0,t.jump.y)[j+1]<=0){sum_j[1,j+1]=0}
      
      if (j_prime!=0){
        for (j_p in 1:j_prime){
          T=tseq[seq]-max(c(0,t.jump.y)[j+1],Ri$stoptime[j_p],na.rm=TRUE)
          Cbj=exp((etam+delta2*vi)*(j_p-1))*(exp(etam+delta2*vi)-1)
          
          if (T <= 0){sum_j[j_p + 1,j + 1]=0}
          if (T > 0){sum_j[j_p + 1,j + 1] = Cb0[j+1] * Cbj * T}
          
        }
      }
      
    }
    sum_2=sum(sum_j)
    S=exp(-exp(etaz*Q + etax%*%Xi+delta1*vi)*sum_2)
    
    S_seq[seq]=S
  }
  
  return(S_seq)
}

