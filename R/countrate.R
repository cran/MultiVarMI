countrate<-function(count.dat) {
  count.rate<-apply(count.dat, 2, function(x) {
    xbar<-mean(x, na.rm=TRUE)
    out<-list(dat=x, rate=xbar)
    return(out)
  })
}

