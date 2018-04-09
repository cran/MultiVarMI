\name{countrate}
\alias{countrate}
\title{Computation of Rates for Count Data}
\description{This function computes the empirical rates for count data.}
\usage{
countrate(count.dat)
}
\arguments{
  \item{count.dat}{A matrix consisting of count variables.}
}
\value{A list of length ncol(count.dat) containing the data and empirical rates for each variable in count.dat.}
\seealso{\code{\link{MI}}, \code{\link{MVN.corr}}}
\examples{
library(PoisBinOrdNonNor)
set.seed(123)
n<-1e4

lambdas<-list(1, 3)

#generate Poisson data
cmat.star <- find.cor.mat.star(cor.mat = .4 * diag(2) + .6, 
                               no.pois = length(lambdas), 
                               pois.list = lambdas)

cntdata <- genPBONN(n, 
                    no.pois = length(lambdas), 
                    cmat.star = cmat.star, 
                    pois.list = lambdas)

#set a sample of the data to missing
cntdata<-apply(cntdata, 2, function(x) {
  x[sample(1:n, size=n/10)]<-NA
  return(x)
})

cntdata<-data.frame(cntdata)
cntinfo<-countrate(cntdata)
}
