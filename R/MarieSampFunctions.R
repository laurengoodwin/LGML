
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


strat.pop.fun<-function(L,N.h,n.h,phat.h,alpha){
  N<-sum(N.h)
  phat.st<-(1/N)*sum(N.h*phat.h)#est of p
  var.phat.st<-sum(((N.h-n.h)/N.h)*((N.h/N)^2)*((phat.h*(1-phat.h))/(n.h-1)))# est of p var
  se.phat.st<-sqrt(var.phat.st)#se of p
  CI<-phat.st+c(-1,1)*qt((1-alpha/2),sum(n.h)-L)*se.phat.st
  cat("phat.st=", phat.st, "SE=", se.phat.st,"CI=",CI, "\n")
}


