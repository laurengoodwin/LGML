
##unbiased_cluster gives an estimate of the population total and a standard error. In addition it will 
# at the given (1-alpha)*100% confidence level, the mean y value per secondary unit, and mean y value per 
# primary unit. The function requires inputs of 
# y.i = the total of the ys in the ith primary unit
# M.i = the number of secondary units in the ith primary unit.
# M = the total number of secondary units in the population
# N = the number of primary units or clusters in the population.
# n = the number of primary units in the sample.


unbiased_cluster<-function(y.i,M.i,M,N,n,alpha){
tau.hat<- (N/n)*sum(y.i)
su.sq<-(1/(n-1))*sum((y.i-mean(y.i))^2)
se.tauhat<- sqrt((N*(N-n))*(su.sq/n))
mu.hat<-tau.hat/M
se.mu.hat<-se.tauhat/M
mu.hat.1<-tau.hat/N
se.mu.hat.1<-tau.hat/N
CI.tau<-tau.hat+c(-1,1)*qt((1-alpha/2),n-1)*se.tauhat
CI.mu<-mu.hat+c(-1,1)*qt((1-alpha/2),n-1)*se.mu.hat
CI.mu1<-mu.hat.1+c(-1,1)*qt((1-alpha/2),n-1)*se.mu.hat.1
cat("tau.hat=", tau.hat, "SE=", se.tauhat,"CI=",CI.tau, "\n")
cat("mu.hat=", mu.hat, "SE=", se.mu.hat,"CI=",CI.mu, "\n")
cat("mu.hat.1=", mu.hat.1, "SE=", se.mu.hat.1,"CI=",CI.mu1, "\n")
}


##ratio_cluster gives an estimate of the population total and a standard error. In addition it will 
# at the given (1-alpha)*100% confidence level, the mean y value per secondary unit, and mean y value per 
# primary unit. The function requires inputs of 
# y.i = the total of the ys in the ith primary unit
# M.i = the number of secondary units in the ith primary unit.
# M = the total number of secondary units in the population
# N = the number of primary units or clusters in the population.
# n = the number of primary units in the sample.

ratio_cluster<-function(y.i,M.i,M,N,n,alpha){
  tau.hatr<- M*(sum(y.i)/sum(M.i))
  mu.hatr<-tau.hatr/M
  M.bar<- (1/n)*sum(M.i)
  sr.sq<-(sum((y.i-mu.hatr*M.i)^2)/(n-1))
  se.tauhatr<-sqrt(N*(N-n)*(sr.sq/n))
  se.muhatr<-se.tauhatr/M
  CI.taur<-tau.hatr+c(-1,1)*qt((1-alpha/2),n-1)*se.tauhatr
  CI.mur<-mu.hatr+c(-1,1)*qt((1-alpha/2),n-1)*se.muhatr
  cat("tau.hat.r=", tau.hatr, "SE=", se.tauhatr,"CI=",CI.taur, "\n")
  cat("mu.hat.r=", mu.hatr, "SE=", se.muhatr,"CI=",CI.mur, "\n")
}

