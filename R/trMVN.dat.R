trMVN.dat<-function(indat, ord.mps=NULL, nct.sum=NULL, count.rate=NULL) {
  if(is.null(ord.mps) & is.null(nct.sum) & is.null(count.rate)) {
    stop('Ord.mps, nct.sum, and count.rate are all set to NULL.')
  }
  if(any(unlist(lapply(ord.mps, sum))>(1+.Machine$double.eps^0.5)) | any(unlist(lapply(ord.mps, sum))<(1-.Machine$double.eps^0.5))) { #tolerance added for use across platforms
    stop('Marginal probabilities in each element of ord.mps must sum to 1.')
  }
  if(any(count.rate<=0)) {
    stop("Rates for count variables must be greater than 0.")
  }
  if(any(nct.sum['Variance',]<=0)) {
    stop("Variance for continuous variables must greater than 0.")
  }
  
  tr.dat<-lapply(indat, function(x) {
    dat.new<-NULL
    
    if(!is.null(nct.sum)) {
      nct.tvars<-colnames(nct.sum)
      for(ivar in nct.tvars) {
        
        #get normalized variable
        ivar.info<-nct.sum[,ivar]
        if(is.null(names(ivar.info))) {
          names(ivar.info)<-rownames(nctsum)
        }
        
        ivar.dat<-x[,ivar]
        ivar.new<-ivar.info['a']+ivar.dat*ivar.info['b']+ivar.dat^2*ivar.info['c']+ivar.dat^3*ivar.info['d']
        
        #mean and sd
        ivar.new<-ivar.new*sqrt(ivar.info['Variance'])+ivar.info['Mean']
        
        #prepare for output
        ivar.new<-matrix(ivar.new)
        colnames(ivar.new)<-ivar
        dat.new<-cbind(dat.new, ivar.new)
      }
    }
    
    if(!is.null(ord.mps)) {
      ord.tvars<-names(ord.mps)
      
      for(ivar in ord.tvars) {
        ivar.mps<-ord.mps[[ivar]]
        ivar.cats<-as.numeric(names(ivar.mps))
        ivar.dat<-x[,ivar]
        ivar.new<-ordY(ivar.mps, ivar.cats, y=ivar.dat)$x
        
        #prepare for output
        ivar.new<-matrix(ivar.new)
        colnames(ivar.new)<-ivar
        dat.new<-cbind(dat.new, ivar.new)
      }
    } 
    
    if(!is.null(count.rate)) {
      count.tvars<-names(count.rate)
      for(ivar in count.tvars) {
        ivar.rate<-count.rate[ivar]
        ivar.dat<-x[,ivar]
        
        U<-pnorm(ivar.dat)
        ivar.new<-qpois(U, ivar.rate)
        
        #prepare for output
        ivar.new<-matrix(ivar.new)
        colnames(ivar.new)<-ivar
        dat.new<-cbind(dat.new, ivar.new)       
      }
    }
    return(dat.new)
  }) 
  return(tr.dat)
}