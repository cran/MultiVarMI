\name{nctsum}
\alias{nctsum}
\title{Summary and Standardization of Continuous Data}
\description{This function calculates mean, variance, skewness, excess kurtosis, and Fleishman coefficients for continuous data and also standardizes each variable.}
\usage{
nctsum(nct.dat)
}
\arguments{
  \item{nct.dat}{A data frame consisting of continuous variables.}
}
\value{A list of length ncol(nct.dat) containing the standardized data and summary statistics for each variable in nct.dat.}
\references{
Fleishman A.I. (1978). A method for simulating non-normal distributions. \emph{Psychometrika}, \bold{43(4)}, 521-532.
}
\seealso{\code{\link{MI}}, \code{\link{MVN.corr}}}
\examples{
library(PoisBinOrdNonNor)
set.seed(123)
n<-1e4

#first four moments for each simulated variable
moms<-list(c(-1, 1, 0, 1), c(0, 3, 0, 2)) 
  
#generate continuous data
cmat.star <- find.cor.mat.star(cor.mat = .8 * diag(2) + .2, 
                               no.nonn = 2, 
                               nonn.list = moms)

nctdata <- genPBONN(n,
                    no.nonn = length(moms),
                    cmat.star = cmat.star, 
                    nonn.list = moms)

#set a sample of each variable to missing
nctdata<-apply(nctdata, 2, function(x) {
  x[sample(1:n, size=n/10)]<-NA
  return(x)
})

nctdata<-data.frame(nctdata)
nctinfo<-nctsum(nctdata)
}
