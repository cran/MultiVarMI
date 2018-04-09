MVN.corr<-function(indat, var.types, ord.mps=NULL, nct.sum=NULL, count.rate=NULL) {
  if(any(var.types=="O") & is.null(ord.mps)) {
    stop('Marginal probabilities must be provided for ordinal variables.')
  }
  if(any(var.types=="C") & is.null(count.rate)) {
    stop('Rates must be provided for count variables.')
  }  
  if(any(var.types=="NCT") & is.null(nct.sum)) {
    stop('Summary must be provided for continuous variables.')
  }
  if(is.null(ord.mps) & is.null(nct.sum) & is.null(count.rate)) {
    stop('Ord.mps, nct.sum, and count.rate are all set to NULL.')
  }
  if(any(unlist(lapply(ord.mps, sum))>(1+.Machine$double.eps^0.5)) | any(unlist(lapply(ord.mps, sum))<(1-.Machine$double.eps^0.5))) { #tolerance added for use across platforms
    stop('Marginal probabilities in each element of ord.mps must sum to 1.')
  }
  if(any(count.rate<=0)) {
    stop("Rates for count variables must be greater than 0.")
  }
  
  #assign variable name to type of variable
  names(var.types)<-colnames(indat)
  
  #create a list with every combination of 2 columns
  get.pairs<-combn(colnames(indat), 2)
  
  list.dfs<-apply(get.pairs, 2, function(x) {
    dat.pair<-indat[,colnames(indat) %in% c(x[1],x[2])]
    type.pair<-paste(sort(var.types[names(var.types) %in% c(x[1],x[2])]), collapse='')
    return(list(type.pair, dat.pair))
  })
  
  all.corr<-lapply(list.dfs, function(x) {
    pair.type<-x[[1]]
    pair.dat<-x[[2]]
    pair<-paste(names(pair.dat), collapse='/')
    pair.corr<-cor(x=pair.dat[,1], y=pair.dat[,2], use="pairwise.complete.obs")
    

    if(pair.type=="OO") {  #correlation for the O-O pairs
      
      pair.mps<-ord.mps[names(ord.mps) %in% colnames(pair.dat)]
      
      out.corr<-ophi2corrZ(ophi=pair.corr, 
                           p1=pair.mps[[1]], 
                           p2=pair.mps[[2]])
      
    } else if(pair.type=="CC") { #correlation for the C-C pairs
      
      c.mat<-matrix(c(1, pair.corr, pair.corr, 1), byrow=TRUE, nrow=2, ncol=2)
      pair.rates<-count.rate[names(count.rate) %in% colnames(pair.dat)]
      
      out.corr<-intercor.PP(lamvec=pair.rates, 
                            cmat=c.mat)[1,2]
      
    } else if(pair.type=="NCTNCT") { #correlation for NCT-NCT
      
      pair.sum<-nct.sum[,colnames(pair.dat)]
      skew.vec<-pair.sum['Skewness',]
      kurto.vec<-pair.sum['Excess Kurtosis',]
      
      out.corr<-corrY2corrZ(corrY=pair.corr, skew.vec=skew.vec, kurto.vec=kurto.vec)

    } else if(pair.type=="NCTO") { #correlation for NCT-O
      
      p<-ord.mps[names(ord.mps) %in% colnames(pair.dat)][[1]]
      cats<-as.numeric(names(p))
      
      X <- rnorm(1e5, 0, 1) 
      Y <- rnorm(1e5, 0, 1) 
      
      XO <- ordY(mp = p, cat = cats, y = X)$x #ordinalize X
      
      cor_ON <-cor(XO[order(XO)], Y[order(Y)]) / cor(X[order(X)], Y[order(Y)])
      cor_NCTN<-pair.corr/cor_ON

      b<-nct.sum['b', colnames(nct.sum) %in% colnames(pair.dat)]
      d<-nct.sum['d', colnames(nct.sum) %in% colnames(pair.dat)]
      
      out.corr<-cor_NCTN/(b+3*d)
      
    } else if(pair.type=="CNCT") { #correlation for NCT-C
      
      C.rate<-count.rate[names(count.rate) %in% colnames(pair.dat)]
      N.fc<-nct.sum[c('a','b','c','d'), colnames(nct.sum) %in% colnames(pair.dat)]
        
      out.corr<-intercor.NNP(lamvec=C.rate,
                             cmat=matrix(pair.corr),
                             pmat=matrix(N.fc, nrow=1))
      
    } else if(pair.type=="CO") { #correlation for O, C
      
      p<-ord.mps[names(ord.mps) %in% colnames(pair.dat)][[1]]
      cats<-as.numeric(names(p))
      
      X <- rnorm(1e5, 0, 1) 
      Y <- rnorm(1e5, 0, 1) 
      
      XO <- ordY(mp = p, cat = cats, y = X)$x #ordinalize X
      
      cor_ON <- cor(XO[order(XO)], Y[order(Y)]) / cor(X[order(X)], Y[order(Y)])
      cor_CN <- pair.corr/cor_ON
      
      rate<-count.rate[names(count.rate) %in% colnames(pair.dat)]
      
      U<-pnorm(Y)
      YC<-qpois(U, rate)
        
      cor_YYC<-cor(Y, YC)
      
      out.corr<-cor_CN/cor_YYC
    }  
    names(out.corr)<-pair
    return(out.corr)
  })
  
  #make correlation matrix
  nvar<-ncol(indat)
  cor.mat <- diag(nvar)
  
  #add colnames/rownames to cor.mat
  colnames(cor.mat)<-rownames(cor.mat)<-colnames(indat)
  
  cor.mat[lower.tri(cor.mat)] <- unlist(all.corr)
  cor.mat[upper.tri(cor.mat)] <- t(cor.mat)[upper.tri(t(cor.mat))]
  
  #check if correlation matrix is positive-definite
  if(!is.positive.definite(cor.mat)) {
    cor.mat<-matrix(nearPD(x=cor.mat, corr=TRUE)$mat, nrow=nvar, ncol=nvar)
  }
  
  return(cor.mat)
}

  

