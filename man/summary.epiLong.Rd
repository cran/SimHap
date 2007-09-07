\name{summary.epiLong}
\alias{summary.epiLong}
\alias{print.summary.epiLong}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing epidemiological analysis models for longitudinal data}
\description{
  Summary method for objects of class \code{epiLong}
}
\usage{
\method{summary}{epiLong}(object, ...)
\method{print}{summary.epiLong}(x, digits = max(3, getOption("digits") - 3), 
	signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{object}{an object of class \code{epiLong}, a result of a call to \code{epi.long}.}
  \item{x}{an object of class \code{summary.epiLong}, the result of a call to \preformatted{summary.epiLong.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.epiLong} returns an object of class \code{summary.epiLong}. Some components taken from \code{summary.lme}. A list with components
  
  \item{call}{the formula call.}
  \item{terms}{terms attribute of the formula called in \code{epi.long}.}
  \item{fixDF}{as per \code{lmeObject}. A list with components \code{X} and \code{terms} specifying the denominator degrees of freedom for, respectively, t-tests for the individual fixed effects and F-tests for the fixed-effects terms in the models.}
  \item{sigma}{as per \code{lmeObject}. The estimated within-group error standard deviation.}
  \item{groups}{as per \code{lmeObject}. A data frame with the grouping factors as columns. The grouping level increases from left to right.}
  \item{dims}{as per \code{lmeObject}. A list with basic dimensions used in the lme fit, including the components \code{N} - the number of observations in the data, \code{Q} - the number of grouping levels, \code{qvec} - the number of random effects at each level from innermost to outermost (last two values are equal to zero and correspond to the fixed effects and the response), \code{ngrps} - the number of groups at each level from innermost to outermost (last two values are one and correspond to the fixed effects and the response), and \code{ncol} - the number of columns in the model matrix for each level of grouping from innermost to outermost (last two values are equal to the number of fixed effects and one).}
  \item{method}{as per \code{lmeObject}. The estimation method: either ``ML" for maximum likelihood, or ``REML" for restricted maximum likelihood.}
  \item{residuals}{as per \code{lmeObject}. If more than five observations are used in the \code{lme} fit, a vector with the minimum, first quartile, median, third quartile, and maximum of the innermost grouping level residuals distribution; else the innermost grouping level residuals.}
  \item{coefficients}{summarized results from fitted model, including coefficients, standard errors and p-values.}
  \item{fixed_formula}{fixed effects formula used in \code{epi.long}.}
  \item{random_formula}{random effects formula used in \code{epi.long}.}
  \item{AIC}{Akaike Information Criterion for the model fitted in \code{epi.long}.}
  \item{corStruct}{the correlation structure specified in \code{epi.long}.}
  \item{modelStruct}{as per \code{lmeObject}. An object inheriting from class \code{lmeStruct}, representing a list of mixed-effects model components, such as \code{reStruct}, \code{corStruct}, and \code{varFunc} objects.}

}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

}
\author{Pamela A. McCaskie}

\seealso{\code{\link{epi.long}}}
\examples{

data(longPheno.dat)
mymodel <- epi.long(fixed=fev1f~height, random=~1|id, pheno=longPheno.dat, 
	form=~year|id)
summary(mymodel)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

