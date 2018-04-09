MVN.dat<-function(ord.info=NULL, nct.info=NULL, count.info=NULL) {
  if(is.null(ord.info) & is.null(nct.info) & is.null(count.info)) {
    stop('Ord.info, nct.info, and count.info are all set to NULL.')
  }
  
  if(!is.null(ord.info)) {
    #transform ordinal variables to normal variables
    ord.mvn<-lapply(ord.info, function(x) {
      tdat<-x$dat
      ndat<-rep(NA, length(tdat))
      
      cps<-unlist(mps2cps(list(x$mps)))
      mvn.cut<-qnorm(cps)
      mvn.samp<-rnorm(5e6)
      
      cnts<-table(x$dat)
      
      for(icat in 1:length(cnts)) {
        if(icat==1) {
          cati<-sample(mvn.samp[which(mvn.samp<=mvn.cut[icat])], size=cnts[icat], replace=TRUE)
        } else if(icat==length(cnts)) {
          cati<-sample(mvn.samp[which(mvn.samp>mvn.cut[(icat-1)])], size=cnts[icat], replace=TRUE)
        } else {
          cati<-sample(mvn.samp[which(mvn.samp>mvn.cut[(icat-1)] & mvn.samp<=mvn.cut[icat])], size=cnts[icat], replace=TRUE)
        }
        ndat[which(tdat==names(cnts[icat]))]<-cati
      }
      
      return(ndat)
    })
    ord.mvn.df<-do.call(cbind, ord.mvn)
  } else {
    ord.mvn.df<-NULL
  }
  
  if(!is.null(nct.info)) {
   nct.mvn<-lapply(nct.info, function(x) {
     tdat<-x$dat.S
     ndat<-numeric()
     for(i in 1:length(tdat)) {
       if(is.na(tdat[i])) {
         ndat[i]<-NA 
       } else {
         sol<-polyroot(c(x$summ['a']-tdat[i], x$summ['b'], x$summ['c'], x$summ['d']))
         if(length(sol)==1) {
           ndat[i]<-Re(sol)[abs(Im(sol)) < 1e-6] 
         } else {
           ndat[i]<-sample(Re(sol)[abs(Im(sol)) < 1e-6], 1) #more than 1 solution, randomly sample
         }
                 
       }
     }
     return(ndat)
   })
   
   nct.mvn.df<-do.call(cbind, nct.mvn)
  } else {
    nct.mvn.df<-NULL
  }   
  
  if(!is.null(count.info)) {
    count.mvn<-lapply(count.info, function(x) {
      tdat<-x$dat
      lambda<-x$rate
      
      #get percentiles of observed counts for given lambda
      U<-ppois(sort(unique(tdat[!is.na(tdat)])), lambda)
      
      #get normal quantiles for count variable percentiles
      mvn.cut<-qnorm(U)
      
      cnts<-table(tdat)
      
      mvn.samp<-rnorm(5e6)
      ndat<-rep(NA, length(tdat))
      
      for(icat in 1:length(cnts)) {
        if(icat==1) {
          cati<-sample(mvn.samp[which(mvn.samp<=mvn.cut[icat])],
                       size=cnts[icat],
                       replace=TRUE)
        } else if(icat==length(cnts)) {
          if(length(mvn.samp[which(mvn.samp>mvn.cut[icat])])>0){
            cati<-sample(mvn.samp[which(mvn.samp>mvn.cut[icat])],
                         size=cnts[icat],
                         replace=TRUE)
          } else {
            #maximum value not in sampled normal distribution
            cati<-sample(max(mvn.cut),
                         size=cnts[icat],
                         replace=TRUE)
          }

          
        } else {
          cati<-sample(mvn.samp[which(mvn.samp>mvn.cut[(icat-1)] & mvn.samp<=mvn.cut[icat])],
                       size=cnts[icat],
                       replace=TRUE)
        }
        ndat[which(tdat==names(cnts[icat]))]<-cati
      }
      return(ndat)
    })
    
    
    count.mvn.df<-do.call(cbind, count.mvn)
  } else {
    count.mvn.df<-NULL
  } 
  
  mvn<-cbind(nct.mvn.df, ord.mvn.df, count.mvn.df)
  return(mvn)
}