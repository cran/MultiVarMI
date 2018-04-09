\name{MI}
\alias{MI}
\title{Bayesian Multiple Imputation for Multivariate Data}
\description{This function implements the multiple imputation framework as described in Demirtas (2017) "A multiple imputation framework for massive multivariate
data of different variable types: A Monte-Carlo technique."}
\usage{
MI(dat, var.types, m)
}
\arguments{
  \item{dat}{A data frame containing multivariate data with missing values.}
  \item{var.types}{The variable type corresponding to each column in dat, taking values of "NCT" for continuous data, "O" for ordinal or binary data, or "C" for count data.}
  \item{m}{The number of stochastic simulations in which the missing values are replaced.}
}
\value{A list containing m imputed data sets.}
\references{
Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. \emph{The American Statistician}, \bold{65(2)}, 104-109.

Demirtas, H., Hedeker, D., and Mermelstein, R. J. (2012). Simulation of massive public health data by power polynomials. \emph{Statistics in Medicine}, \bold{31(27)}, 3337-3346.

Demirtas, H. (2016). A note on the relationship between the phi coefficient and the tetrachoric correlation under nonnormal underlying distributions. \emph{The American Statistician}, \bold{70(2)}, 143-148.

Demirtas, H. and Hedeker, D. (2016). Computing the point-biserial correlation under any underlying continuous distribution. \emph{Communications in Statistics-Simulation and Computation}, \bold{45(8)}, 2744-2751.

Demirtas, H., Ahmadian, R., Atis, S., Can, F.E., and Ercan, I. (2016). A nonnormal look at polychoric correlations: modeling the change in correlations before and after discretization. \emph{Computational Statistics}, \bold{31(4)}, 1385-1401.

Demirtas, H. (2017). A multiple imputation framework for massive multivariate data of different variable types: A Monte-Carlo technique. \emph{Monte-Carlo Simulation-Based Statistical Modeling}, edited by Ding-Geng (Din) Chen and John Dean Chen, Springer, 143-162.

Ferrari, P.A. and Barbiero, A. (2012). Simulating ordinal data. \emph{Multivariate Behavioral Research}, \bold{47(4)}, 566-589.

Fleishman A.I. (1978). A method for simulating non-normal distributions. \emph{Psychometrika}, \bold{43(4)}, 521-532.

Vale, C.D. and Maurelli, V.A. (1983). Simulating multivariate nonnormal distributions. \emph{Psychometrika}, \bold{48(3)}, 465-471.
}
\seealso{\code{\link{MVN.corr}}, \code{\link{MVN.dat}}, \code{\link{trMVN.dat}}}
\examples{
\donttest{
library(PoisBinOrdNonNor)
set.seed(1234)
n<-1e5
lambdas<-list(1, 3) #2 count variables
mps<-list(c(.2, .8), c(.6, 0, .3, .1)) #1 binary variable, 1 ordinal variable with skip pattern
moms<-list(c(-1, 1, 0, 1), c(0, 3, 0, 2)) #2 continuous variables

################################################
#Generate Poisson, Ordinal, and Continuous Data#
################################################
#get intermediate correlation matrix
cmat.star <- find.cor.mat.star(cor.mat = .8 * diag(6) + .2, #all pairwise correlations set to 0.2
                               no.pois = length(lambdas),
                               no.ord = length(mps),
                               no.nonn = length(moms),
                               pois.list = lambdas,
                               ord.list = mps,
                               nonn.list = moms)

#generate dataset
mydata <- genPBONN(n,
                   no.pois = length(lambdas),
                   no.ord = length(mps),
                   no.nonn = length(moms),
                   cmat.star = cmat.star,
                   pois.list = lambdas,
                   ord.list = mps,
                   nonn.list = moms)

cor(mydata)
apply(mydata, 2, mean)

#Make 10 percent of each variable missing completely at random
mydata<-apply(mydata, 2, function(x) {
  x[sample(1:n, size=n*0.1)]<-NA
  return(x) 
  }
)


#Create 5 imputed datasets
mydata<-data.frame(mydata)
mymidata<-MI(dat=mydata,
             var.types=c('C', 'C', 'O', 'O', 'NCT', 'NCT'),
             m=5)

#get the means of each variable for the m imputed datasets  
do.call(rbind, lapply(mymidata, function(x) apply(x, 2, mean)))

#get m correlation matrices of for the m imputed dataset
lapply(mymidata, function(x) cor(x))

#Look at the second imputed dataset
head(mymidata$dataset2)

##run a linear model on each dataset and extract coefficients
mycoef<-lapply(mymidata, function(x) {
  fit<-lm(X6~., data=data.frame(x))
  fit.coef<-coef(fit)
  return(fit.coef)
})

do.call(rbind, mycoef)
}
}
