\name{MVN.dat}
\alias{MVN.dat}
\title{Computation of Normal Scores for Multivariate Data}
\description{This function assigns a normal score to binary and ordinal variables using normal quantiles in this appropriate range dictated by marginal proportions; a normal score to count variables based on the equivalence of CDFs of Poisson and normal distribution in the appropriate range dictated by the rate parameters; and a normal score for each continuous measurement by finding the normal root in the Fleishman equation.}
\usage{
MVN.dat(ord.info=NULL, nct.info=NULL, count.info=NULL)
}
\arguments{
  \item{ord.info}{A list containing binary and ordinal data and corresponding marginal probabilities as packaged in \code{ordmps}. Default is NULL.}
  \item{nct.info}{A list containing standardized continuous data and corresponding summary statistics for continuous variables as packaged in \code{nctsum}. Default is NULL.}
  \item{count.info}{A list containing count data and corresponding rates as packaged in \code{countrate}. Default is NULL.}
}
\value{A matrix containing normal scores for each variable input.}
\references{
Fleishman A.I. (1978). A method for simulating non-normal distributions. \emph{Psychometrika}, \bold{43(4)}, 521-532.
}
\seealso{\code{\link{MI}}, \code{\link{ordmps}}, \code{\link{nctsum}}, \code{\link{countrate}}}
\examples{
library(PoisBinOrdNonNor)
n<-1e4
lambdas<-list(1)
mps<-list(c(.2, .8))
moms<-list(c(-1, 1, 0, 1))
  
#generate Poisson, ordinal, and continuous data
cmat.star <- find.cor.mat.star(cor.mat = .8 * diag(3) + .2, 
                               no.pois = length(lambdas), 
                               no.ord = length(mps),
                               no.nonn = length(moms), 
                               pois.list = lambdas, 
                               ord.list = mps, 
                               nonn.list = moms)

mydata <- genPBONN(n, 
                   no.pois = length(lambdas), 
                   no.ord = length(mps), 
                   no.nonn = length(moms),
                   cmat.star = cmat.star, 
                   pois.list = lambdas,
                   ord.list = mps, 
                   nonn.list = moms)

#set a sample of each variable to missing
mydata<-apply(mydata, 2, function(x) {
  x[sample(1:n, size=n/10)]<-NA
  return(x)
})

mydata<-data.frame(mydata)

#get information for use in function
count.info<-countrate(count.dat=data.frame(mydata[,c('X1')]))
ord.info<-ordmps(ord.dat=data.frame(mydata[,c('X2')]))
nct.info<-nctsum(nct.dat=data.frame(mydata[,c('X3')]))

mvn.dat<-MVN.dat(ord.info=ord.info,
                 nct.info=nct.info,
                 count.info=count.info) #outputs in order of continuous, ordinal, count
                   
}
