nctsum<-function(nct.dat) {
  
  #get mean, var, skew, excess kurtosis
  nct.sum<-apply(nct.dat, 2, function(x) {
    xbar<-mean(x, na.rm=TRUE)
    s2<-var(x, na.rm=TRUE)
    skew<-skewness(x, na.rm=TRUE)
    exkurt<-kurtosis(x, na.rm=TRUE)-3
    coef.mat <- Fleishman.coef.NN(skew.vec = skew, kurto.vec = exkurt)
    summ<-c(xbar, s2, skew, exkurt, coef.mat)
    names(summ)<-c('Mean', 'Variance', 'Skewness', 'Excess Kurtosis', 'a', 'b', 'c', 'd')
    
    #standardize continuous data
    x.S<-(x-xbar)/(sd(x, na.rm=TRUE))
    
    out<-list(dat.S=x.S, summ=summ)
    return(out)
  })
  
  return(nct.sum)
}