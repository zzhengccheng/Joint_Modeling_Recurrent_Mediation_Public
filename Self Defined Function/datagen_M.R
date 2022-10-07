#Function to generate M0 and M1
datagen_M<-function(X,vi,betaz,betax,lambda.m,t.jump.m,tau){
  
  ####Simulate M^z(t)
  ID=1:n
  beta.all=c(betaz,betax,1)
  risk0=c(exp(cbind(0,X,vi)%*%beta.all))
  risk1=c(exp(cbind(1,X,vi)%*%beta.all))
  ####generate survival time until all individual has follow-up at lest time tau
  M0=mysimRec(t.jump.m,lambda.m,tau,risk0,recurrent=TRUE)
  M1=mysimRec(t.jump.m,lambda.m,tau,risk1,recurrent=TRUE)
  M0$M=0
  M1$M=1
  rbind(M0,M1)
}