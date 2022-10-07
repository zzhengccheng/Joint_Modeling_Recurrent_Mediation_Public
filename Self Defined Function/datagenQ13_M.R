####Function used to generate recurrent events M_Q1 and M_Q3
datagenQ13_M<-function(X,Q1,Q3,vi,betaz,betax,tau,lambda.m,t.jump.m){
  
  ####Simulate M^z(t)
  ID=1:n
  beta.all=c(betaz,betax,1)
  risk_Q1=c(exp(cbind(Q1,X,vi)%*%beta.all))
  risk_Q3=c(exp(cbind(Q3,X,vi)%*%beta.all))
  ####generate survival time until all individual has follow-up at lest time tau
  M_Q1=mysimRec(t.jump.m,lambda.m,tau,risk_Q1,recurrent=TRUE)
  M_Q3=mysimRec(t.jump.m,lambda.m,tau,risk_Q3,recurrent=TRUE)
  M_Q1$M=Q1
  M_Q3$M=Q3
  rbind(M_Q1,M_Q3)
}