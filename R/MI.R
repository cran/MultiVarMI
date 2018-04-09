MI<-function(dat, var.types, m) {
  if(ncol(dat)!=length(var.types)) {
    stop("Variable type must be specified for every column.") 
  }
  if(!any(var.types %in% c('NCT', 'O', 'C'))) {
    stop("Variable types include NCT, O, and C.")
  }
  if(m<=0 | m %% 1 !=0) {
    stop('Number of iterations m must be a positive integer.')
  }
  if(length(dat[is.na(dat)])==0) {
    stop("There are no missing variables in dat.")
  }

  #assign variable name to type of variable
  names(var.types)<-colnames(dat)
 
  #get attributes for each column based on observed data
  ord.dat<-dat[which(var.types=='O')]
  count.dat<-dat[which(var.types=='C')]
  nct.dat<-dat[which(var.types %in% c('NCT'))] 
  
  if(ncol(ord.dat)>0) { #get observed probabilities
    ord.info<-ordmps(ord.dat=ord.dat)
    mps<-lapply(ord.info, "[[", 2) #list
  } else {
    ord.info<-NULL
    mps<-NULL
  }
  
  if(ncol(nct.dat)>0) {
    nct.info<-nctsum(nct.dat=nct.dat)

    #replace continuous data with standardized form
    nct.dat.S<-sapply(nct.info, "[[", 1) #get the columns of standardized continuous data
    
    #remove NCT columns in original dataset and replace with standardized form
    dat[,colnames(nct.dat.S)]<-nct.dat.S[,colnames(nct.dat.S)]
    
    #store moments (mean, variance, skewness, excess kurtosis) and Fleishman coefficients
    allnctsum<-sapply(nct.info, "[[", 2) #matrix
  } else {
    nct.info<-NULL
    allnctsum<-NULL
  }
  
  if(ncol(count.dat)>0) {
    count.info<-countrate(count.dat=count.dat)
    rates<-sapply(count.info, "[[", 2) #numeric vector
  } else {
    count.info<-NULL
    rates<-NULL
  }
  
  #get NORMAL correlation of each pair (list.dfs[[i]][[1]] gives type.pair, list.dfs[[i]][[2]] gives data)
  mvn.cmat<-MVN.corr(indat=dat,
                     var.types=var.types,
                     ord.mps=mps,
                     nct.sum=allnctsum,
                     count.rate=rates)
  
  #Create MVN data
  mvn.dat<-MVN.dat(ord.info=ord.info,
                   nct.info=nct.info,
                   count.info=count.info)
 
  #Bayesian MI creating m datasets
  mvn.means<-rep(0, ncol(mvn.dat))
  mvn.sd<-rep(1, ncol(mvn.dat))
  
  rngseed(1234567)
  s<-prelim.norm(mvn.dat)
  thetahat<-makeparam.norm(s, thetalist=list(mvn.means, mvn.sd, mvn.cmat))
  
  mdat<-replicate(m, mvn.dat, simplify = FALSE)
  
  mi.dat<-lapply(mdat, function(x) {
    da.dat<-da.norm(s, start=thetahat, steps=100, showits=FALSE, return.ymis=TRUE)
    x[is.na(x)]<-da.dat$ymis
    return(x)
  })
  
  #Backtransform m datasets
  final.dat<-trMVN.dat(indat=mi.dat, 
                       ord.mps=mps,
                       nct.sum=allnctsum,
                       count.rate=rates)
  
  final.dat<-lapply(final.dat, function(x) {
    return(x[,colnames(dat)])
  })
  
  names(final.dat)<-paste0("dataset", 1:m)
  
  return(final.dat)
}










