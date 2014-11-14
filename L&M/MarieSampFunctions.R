#Q2 from Midterm
n.plants<-n.plants<-c(0,22,1,12,4,21,77,27,4,23,19,63,55,39,17,7,18,3,12,5)#y
area<-c(22,27,33,38.4,41,45,34,34,41,31,37,42,30.5,30,27,3,13,12,14,13)#x
N<-157
#a
SRS.est<-function(y,alpha,N=NA){
  # Estimate mean and total.
  # y is the data
  # N is the population size (assumed infinite if no value is given)
  n<-length(y)
  fpc<-1
  if(!is.na(N)) fpc<-(N-n)/N
  ybar<-mean(y)
  SE.ybar<-sqrt(fpc * var(y)/n)
  if(length(y) < 30) warning("sample size too small-CI may not be valid")
  y.bar.CI<-ybar+c(-1,1)*qt((1-alpha/2),n-1)*SE.ybar
  cat("ybar=", ybar, "SE=", SE.ybar,"CI=",y.bar.CI, "\n")
  if(!is.na(N))
    cat("tau.hat=", N*ybar, "SE=", N*SE.ybar,"CI=",N*ybar+c(-1,1)*qt((1-alpha/2),n-1)*N*SE.ybar, "\n")#est of pop total and se from sample mean
}
SRS.est(n.plants,.05,N)


#b
mux<-5897/157
ratio.est <- function(x, y, mux = NA, N = NA,alpha) {
  # estimate of a ratio and ratio estimate of population mean and total.
  # x is auxiliary variable, y is response, mux is population mean
  # of x (xbar is used if no value is given),
  # N is population size (assumed infinite if no value given),
  if(length(x) != length(y)) stop("x and y must be same length")
  n <- length(x)
  fpc <- 1
  if(!is.na(N))
    fpc <- (N - n)/N
  r <- sum(y)/sum(x)
  sr2 <- (1/(n - 1)) * sum((y - r * x)^2)
  if(is.na(mux)) mx <- mean(x) else mx <- mux
  SE.r<-sqrt((fpc * sr2)/(mx^2 * n))
  if(length(y) < 30) warning("sample size too small-CI may not be valid")
  r.CI<-r+c(-1,1)*qt((1-alpha/2),n-1)*SE.r
  cat("r=", r, " SE=", SE.r,"r.CI=" ,r.CI ,"\n")
  if(!is.na(mux))
    SE.mu<-sqrt((fpc * sr2)/n)
  mu.hat<- r * mux
  mu.CI<- mu.hat+c(-1,1)*qt((1-alpha/2),n-1)*SE.mu
    cat("mu.hat=", mu.hat, " SE=",SE.mu,"mu.CI=",mu.CI , "\n")
  if(!is.na(N) & !is.na(mux))
    tau.hat<-N*r * mux
  SE.tau<-N*sqrt((fpc * sr2)/n)
  tau.CI<- tau.hat+c(-1,1)*qt((1-alpha/2),n-1)*SE.tau
    cat("tau.hat=",tau.hat , " SE.tau=",SE.tau,"tau.CI=",tau.CI , "\n")
}
ratio.est(area,n.plants,mux,N,.05)
================================================================================================================
#Q one from HW5
## Enter ages
y<-c(125,119,83,85,99,117,69,133,154,168,61,80,114,147,122,106,82,88,97,99)
## Enter diameters
x<-c(120,114,79,90,105,79,73,102,117,113,57,80,103,120,92,85,70,107,93,82)/10
regr.est <- function(x, y, mux, N = NA,alpha) {
  # regression estimator of a population mean and total.
  # x is auxiliary variable, y is response, mux is population mean
  # of x, N is population size (assumed infinite if no value given).
  if(length(x) != length(y)) stop("x and y must be same length")
  n <- length(x)
  fpc <- 1
  if(!is.na(N))
    fpc <- (N - n)/N
  ab <- lsfit(x, y)
  mu.hat<-ab$coef[1] + ab$coef[2] * mux
  SE.mu<-sqrt((fpc *sum(ab$residual^2))/(n * (n -2)))
  mu.CI<- mu.hat+c(-1,1)*qt((1-alpha/2),n-1)*SE.mu
  cat("mu.hat=", mu.hat, " SE.mu=", SE.mu,"mu.CI=",mu.CI, "\n")
  if(!is.na(N))
  tau.hat<-N*(ab$coef[1] + ab$coef[2] * mux)
  SE.tau<-N*sqrt((fpc *sum(ab$residual^2))/(n * (n -2)))
  tau.CI<- tau.hat+c(-1,1)*qt((1-alpha/2),n-1)*SE.tau
    cat("tau.hat=",tau.hat, " SE.tau=",SE.tau ,"tau.CI=",tau.CI, "\n")
}
regr.est(x,y,10.3,1132,.05)


===============================================================================================================
  #Strata Functions
  #I checked these with examples for his notes.  See strat PDF pages 1-7
  #Function for mean
  N<-41
L<-3
N.h<-c(20,9,12)
n.h<-c(5,3,4)
ybar.h<-c(1.6,2.8,0.6)
s.2.h<-c(3.3,4.0,2.2)

strat.mean.fun<-function(N,L,N.h,n.h,y.bar.h,s.2.h,alpha){
  mu.hat.st<-(1/N)*sum(N.h*ybar.h)# est of strat mean
  var.hat.mu.hat.st<-sum(((N.h-n.h)/N.h)*((N.h/N)^2)*(s.2.h/n.h))#est of strat var
  se.mu.hat.st<-sqrt(var.hat.mu.hat.st)#se of strat
  CI<-mu.hat.st+c(-1,1)*qt((1-alpha/2),sum(n.h)-L)*se.mu.hat.st
  cat("mu.hat.str=", mu.hat.st, "SE=", se.mu.hat.st,"CI=",CI, "\n")
  tau.hat.st<-sum(N.h*ybar.h)#total est
  var.tau.hat.st<-(N^2)*var.hat.mu.hat.st#est of total var
  se.tau.hat.st<-sqrt(var.tau.hat.st)#se of total
  CI<-tau.hat.st+c(-1,1)*qt((1-alpha/2),sum(n.h)-L)*se.tau.hat.st
  cat("tau.hat.str=", tau.hat.st, "SE=", se.tau.hat.st,"CI=",CI, "\n")
}
strat.mean.fun(N,L,N.h,n.h,y.bar.h,s.2.h,.05)

===============================================================================================================
  #Function for proportion
  N.h<-c(9100,1950,5500,10850,2100,5500,9000)
n.h<-c(636,451,481,611,493,575,588)
phat.h<-c(.38,.27,.18,.19,.36,.13,.26)
L<-7

strat.pop.fun<-function(L,N.h,n.h,phat.h,alpha){
  N<-sum(N.h)
  phat.st<-(1/N)*sum(N.h*phat.h)#est of p
  var.phat.st<-sum(((N.h-n.h)/N.h)*((N.h/N)^2)*((phat.h*(1-phat.h))/(n.h-1)))# est of p var
  se.phat.st<-sqrt(var.phat.st)#se of p
  CI<-phat.st+c(-1,1)*qt((1-alpha/2),sum(n.h)-L)*se.phat.st
  cat("phat.st=", phat.st, "SE=", se.phat.st,"CI=",CI, "\n")
}
strat.pop.fun(L,N.h,n.h,phat.h,.05)

