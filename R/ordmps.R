ordmps<-function(ord.dat) {
  ord.mps<-lapply(ord.dat, function(x) {
    mps<-as.data.frame(table(x)/sum(table(x)))
    ord.p<-mps$Freq
    names(ord.p)<-as.character(mps$x)
    out<-list(dat=x, mps=ord.p)
    return(out)
  })
  return(ord.mps)
}