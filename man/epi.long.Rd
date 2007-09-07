\name{epi.long}
\alias{epi.long}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Epidemiological analysis for longitudinal data}
\description{
  \code{epi.long} is used to fit linear mixed effects models to epidemiological data with longitudinal outcomes.
}
\usage{
epi.long(fixed, random, pheno, cor="corCAR1", value = 0.2, 
	form=~1, sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{as per \code{lme}. A two-sided linear formula object describing the fixed-effects part of the model including SNP parameters, with the response on the left of a \code{~} operator and the terms, separated by \code{+} operators.}
  \item{random}{as per \code{lme}. A one-sided formula of the form \code{~x1+...+xn | g1/.../gm}, with \code{x1+...+xn} specifying the model for the random effects and \code{g1/.../gm} the grouping structure (\code{m} may be equal to 1, in which case no \code{/} is required). The random effects formula will be repeated for all levels of grouping, in the case of multiple levels of grouping.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{cor}{a \code{corStruct} object describing the within-group correlation structure. Available correlation structures are \code{corAR1}, \code{corCAR1}, and \code{corCompSymm}. See the documentation of \code{corClasses} for a description of these. Defaults to \code{corCAR1}.}
  \item{value}{for \code{corAR1} - the value of the lag 1 autocorrelation, which must be between -1 and 1. For \code{corCAR1} - the correlation between two observations one unit of time apart. Must be between 0 and 1. For \code{corCompSymm} - the correlation between any two correlated observations. Defaults to \code{0.2}.}
  \item{form}{a one sided formula of the form \code{~ t}, or \code{~ t | g}, specifying a time covariate \code{t} and,  optionally, a grouping factor \code{g}. A covariate for this correlation structure must be integer valued. When a grouping factor is present in \code{form}, the correlation structure is assumed to apply only to observations within the same grouping level; observations with different grouping levels are assumed to be uncorrelated. Defaults to \code{~ 1}, which corresponds to using the order of the observations in the data as a covariate, and no groups.}
  \item{sub}{an expression representing a subset of the data on which to perform the models.}
}
\details{
\code{cor} will always default to \code{corCAR1} and \code{value} will always default to \code{0.2}. Be sure to change both parameters accordingly if desired. See \code{\link{corClasses}} for more details.
}
\value{
  \code{epi.long} returns an object of class \code{epiLong}.

     The \code{summary} function can be used to obtain and print a
     summary of the results.  

     An object of class \code{epiLong} is a list containing the
     following components:
     
  \item{results}{a table containing the coefficients, standard errors and p-values of the parameter estimates.}
  \item{fixed_formula}{fixed effects formula.}
  \item{random_formula}{random effects formula.}
  \item{fit.lme}{a \code{lme} object fit using formula.}
  \item{ANOD}{analysis of deviance table for the fitted model.}
  \item{logLik}{the log-likelihood for the fitted model.}
  \item{AIC}{Akaike Information Criterion for the model fit using formula.}
  \item{corStruct}{correlation structure used in the fitted model.}
}
\references{
Bates, D.M., Pinheiro, J.C. (1998) Computational methods for multilevel models. Available in PostScript or PDF formats at     http://franz.stat.wisc.edu/pub/NLME/

Box, G.E.P., Jenkins, G.M., Reinsel, G.C. (1994) \emph{Time Series Analysis: Forecasting and Control}, 3rd Edition, Holden-Day.

Davidian, M., Giltinan, D.M. (1995) \emph{Nonlinear Mixed Effects Models for Repeated Measurement Data}, Chapman and Hall.

Laird, N.M., Ware, J.H. (1982) Random-Effects Models for Longitudinal Data, \emph{Biometrics}, 38, 963-974.

Lindstrom, M.J., Bates, D.M. (1988) Newton-Raphson and EM Algorithms for Linear Mixed-Effects Models for Repeated-Measures Data, \emph{Journal of the American Statistical Association}, 83, 1014-1022.

Littel, R.C., Milliken, G.A., Stroup, W.W., Wolfinger, R.D. (1996) SAS Systems for Mixed Models, SAS Institute.

McCaskie, P.A., Carter, K.W, Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

Pinheiro, J.C., Bates, D.M.  (1996) Unconstrained Parametrizations for Variance-Covariance Matrices, \emph{Statistics and Computing}, 6, 289-296.

Pinheiro, J.C., Bates, D.M. (2000) \emph{Mixed-Effects Models in S and S-PLUS}, Springer.

}
\author{Pamela A McCaskie}
\seealso{\code{\link{snp.long}}, \code{\link{haplo.long}}, \code{\link{corClasses}}}
\examples{

data(longPheno.dat)
mymodel <- epi.long(fixed=fev1f~height+weight, random=~1|id, 
	pheno=longPheno.dat, form=~year|id)
summary(mymodel)

# example with a subsetting variable, looking at males only
mymodel <- epi.long(fixed=fev1f~height+weight, random=~1|id, 
	pheno=longPheno.dat, form=~year|id, sub=expression(sex==1))
summary(mymodel)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

