\name{MVN.corr}
\alias{MVN.corr}
\title{Calculation of Intermediate Correlation Matrix}
\description{This function calculates an intermediate correlation matrix for Poisson, ordinal, and continuous random variables, with specified target correlations and marginal properties.}
\usage{
MVN.corr(indat, var.types, ord.mps=NULL, nct.sum=NULL, count.rate=NULL)
}
\arguments{
  \item{indat}{A data frame containing multivariate data. Continuous variables should be standardized.}
  \item{var.types}{The variable type corresponding to each column in dat, taking values of "NCT" for continuous data, "O" for ordinal or binary data, or "C" for count data.}
  \item{ord.mps}{A list containing marginal probabilities for binary and ordinal variables as packaged from output in \code{ordmps}. Default is NULL.}
  \item{nct.sum}{A matrix containing summary statistics for continuous variables as packaged from output in \code{nctsum}. Default is NULL.}
  \item{count.rate}{A vector containing rates for count variables as packaged from output in \code{countrate}. Default is NULL.}
}
\value{The intermediate correlation matrix.}
\references{
Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. \emph{The American Statistician}, \bold{65(2)}, 104-109.

Demirtas, H., Hedeker, D., and Mermelstein, R. J. (2012). Simulation of massive public health data by power polynomials. \emph{Statistics in Medicine}, \bold{31(27)}, 3337-3346.

Demirtas, H. (2016). A note on the relationship between the phi coefficient and the tetrachoric correlation under nonnormal underlying distributions. \emph{The American Statistician}, \bold{70(2)}, 143-148.

Demirtas, H. and Hedeker, D. (2016). Computing the point-biserial correlation under any underlying continuous distribution. \emph{Communications in Statistics-Simulation and Computation}, \bold{45(8)}, 2744-2751.

Demirtas, H., Ahmadian, R., Atis, S., Can, F.E., and Ercan, I. (2016). A nonnormal look at polychoric correlations: modeling the change in correlations before and after discretization. \emph{Computational Statistics}, \bold{31(4)}, 1385-1401.

Ferrari, P.A. and Barbiero, A. (2012). Simulating ordinal data. \emph{Multivariate Behavioral Research}, \bold{47(4)}, 566-589.

Fleishman A.I. (1978). A method for simulating non-normal distributions. \emph{Psychometrika}, \bold{43(4)}, 521-532.

Vale, C.D. and Maurelli, V.A. (1983). Simulating multivariate nonnormal distributions. \emph{Psychometrika}, \bold{48(3)}, 465-471.
}
\seealso{\code{\link{MI}}, \code{\link{MVN.dat}}, \code{\link{ordmps}}, \code{\link{nctsum}}, \code{\link{countrate}}}
\examples{
\donttest{
library(PoisBinOrdNonNor)
n<-1e4
lambdas<-list(1, 3)
mps<-list(c(.2, .8), c(.6, 0, .3, .1))
moms<-list(c(-1, 1, 0, 1), c(0, 3, 0, 2))
  
#generate Poisson, ordinal, and continuous data
cmat.star <- find.cor.mat.star(cor.mat = .8 * diag(6) + .2, 
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
ord.info<-ordmps(ord.dat=mydata[,c('X3', 'X4')])
nct.info<-nctsum(nct.dat=mydata[,c('X5', 'X6')])
count.info<-countrate(count.dat=mydata[,c('X1', 'X2')])

#extract marginal probabilites, continuous properties, and count rates
mps<-sapply(ord.info, "[[", 2)
nctsum<-sapply(nct.info, "[[", 2)
rates<-sapply(count.info, "[[", 2)

#replace continuous with standardized forms
mydata[,c('X5', 'X6')]<-sapply(nct.info, "[[", 1)[,c('X5', 'X6')] 

var.types<-c('C', 'C', 'O', 'O', 'NCT', 'NCT')
mvn.cmat<-MVN.corr(indat=mydata,
                   var.types=var.types,
                   ord.mps=mps,
                   nct.sum=nctsum,
                   count.rate=rates)
}
}
