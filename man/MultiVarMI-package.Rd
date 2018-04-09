\name{MultiVarMI-package}
\alias{MultiVarMI-package}
\alias{MultiVarMI}
\docType{package}
\title{
Multiple Imputation for Multivariate Data
}
\description{
This package implements a Bayesian multiple imputation framework for multivariate data. Most incomplete data sets constist of interdependent binary, ordinal, count, and continuous data. Furthermore, planned missing data designs have been developed to reduce respondent burden and lower the cost associated with data collection. The unified, general-purpose multiple imputation framework described in Demirtas (2017) can be utilized in developing power analysis guidelines for intensive multivariate data sets that are collected via increasingly popular real-time data capture (RTDC) approaches. This framework can accommodate all four major types of variables with a minimal set of assumptions. The data are prepared for multivariate normal multiple imputation for use in the \code{norm} package and subsequently backtransformed to the original distribution.

This package consists of one main function and six auxiliary functions. Multiple imputation can be performed using the function \code{\link{MI}}. While the auxiliary functions are utilized in \code{\link{MI}}, they can be used as stand-alone functions. \code{\link{nctsum}} outputs a list with summary statistics and Fleishman coefficients and standardized forms of each variable. \code{\link{ordmps}} is utilized for ordinal variables and outputs a list with empirical marginal probabilities and the associated observations for each ordinal variable. \code{\link{countrate}} is designed for variables and outputs a list with empirical rates and the associated observations for each ordinal variable. \code{\link{MVN.corr}} calculates the intermediate correlation matrix, \code{\link{MVN.dat}} transforms variables to a standard normal variable, and \code{\link{trMVN.dat}} transforms standard normal variables to ordinal, count, and/or non-normal continuous variables through specified parameters.
}
\details{
\tabular{ll}{
Package: \tab MultiVarMI\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2018-04-08\cr
License: \tab GPL-2 | GPL-3\cr
}}
\author{
Rawan Allozi, Hakan Demirtas

Maintainer: Rawan Allozi <ralloz2@uic.edu>
}
\references{
Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. \emph{The American Statistician}, \bold{65(2)}, 104-109.

Demirtas, H., Hedeker, D., and Mermelstein, R. J. (2012). Simulation of massive public health data by power polynomials. \emph{Statistics in Medicine}, \bold{31(27)}, 3337-3346.

Demirtas, H. (2016). A note on the relationship between the phi coefficient and the tetrachoric correlation under nonnormal underlying distributions. \emph{The American Statistician}, \bold{70(2)}, 143-148.

Demirtas, H. and Hedeker, D. (2016). Computing the point-biserial correlation under any underlying continuous distribution. \emph{Communications in Statistics-Simulation and Computation}, \bold{45(8)}, 2744-2751.

Demirtas, H., Ahmadian, R., Atis, S., Can, F.E., and Ercan, I. (2016). A nonnormal look at polychoric correlations: modeling the change in correlations before and after discretization. \emph{Computational Statistics}, \bold{31(4)}, 1385-1401.

Demirtas, H. (2017). A multiple imputation framework for massive multivariate data of different variable types: A Monte-Carlo technique. \emph{Monte-Carlo Simulation-Based Statistical Modeling}, edited by Ding-Geng (Din) Chen and John Dean Chen, Springer, 143-162.

Ferrari, P.A. and Barbiero, A. (2012). Simulating ordinal data. \emph{Multivariate Behavioral Research}, \bold{47(4)}, 566-589.

Fleishman A.I. (1978). A method for simulating non-normal distributions. \emph{Psychometrika}, \bold{43(4)}, 521-532.

Vale, C.D. and Maurelli, V.A. (1983). Simulating multivariate nonnormal distributions. 
\emph{Psychometrika}, \bold{48(3)}, 465-471.
}
