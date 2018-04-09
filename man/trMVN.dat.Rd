\name{trMVN.dat}
\alias{trMVN.dat}
\title{Transformation of Normal Scores}
\description{This function backtransforms normal scores for ordinal variables using the thresholds determined by the marginal proportions using quantiles of the normal distribution; normal scores for continuous variables by the sum of linear combinations of standard normals using the corresponding Fleishman coefficients; and normal scores for count variables by the inverse cdf matching procedure.}
\usage{
trMVN.dat(indat, ord.mps=NULL, nct.sum=NULL, count.rate=NULL)
}
\arguments{
  \item{indat}{A list of data frames of normal scores to be backtransformed.}
  \item{ord.mps}{A list containing marginal probabilities for binary and ordinal variables as packaged from output in \code{ordmps}. Default is NULL.}
  \item{nct.sum}{A matrix containing summary statistics for continuous variables as packaged from output in \code{nctsum}. Default is NULL.}
  \item{count.rate}{A vector containing rates for count variables as packaged from output in \code{countrate}. Default is NULL.}
}
\value{A list containing backtransformed data.}
\references{
Fleishman A.I. (1978). A method for simulating non-normal distributions. \emph{Psychometrika}, \bold{43(4)}, 521-532.
}
\seealso{\code{\link{MI}}, \code{\link{ordmps}}, \code{\link{nctsum}}, \code{\link{countrate}}}
\examples{
sndat<-data.frame(matrix(rnorm(1e4), ncol=5, nrow=1e4/5))

#ordinal marginal probabilities
m1<-c(0.4, 0.6)
names(m1)<-c(0,1)
m2<-c(0.2, 0.3, 0.5)
names(m2)<-c(0,2,3)
mps<-list(X1=m1, X2=m2)

#count rates
rates<-c(2, 3)
names(rates)<-c('X3', 'X4')

#continuous
nctsum<-data.frame(X5=c(1, 1, -0.31375, 0.82632, 0.31375, 0.02271)) #Weibull(1,1)
rownames(nctsum)<-c('Mean', 'Variance', 'a', 'b', 'c', 'd')

trdat<-trMVN.dat(indat=list(sndat), ord.mps=mps, nct.sum=nctsum, count.rate=rates)
}
