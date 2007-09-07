\name{summary.epiQuant}
\alias{summary.epiQuant}
\alias{print.summary.epiQuant}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing epidemiological analysis models for quantitative outcomes}
\description{
 Summary method for objects of class \code{epiQuant}
}
\usage{
\method{summary}{epiQuant}(object, ...)
\method{print}{summary.epiQuant}(x, digits = max(3, getOption("digits") - 3), 
	signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{epiQuant}, a result of a call to \code{epi.quant}.}
  \item{x}{an object of class \code{summary.epiQuant}, the result of a call to \preformatted{summary.epiQuant.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  \code{summary.epiQuant} returns an object of class summary.epiQuant, a list with components
  
  \item{call}{the formula call.}
  \item{terms}{terms attribute of the formula called in \code{epi.quant}.}
  \item{na.action}{method used for missing data.}
  \item{coefficients}{summarized results from fitted model, including coefficients, standard errors and p-values.}
  \item{df}{residual degrees of freedom.}
  \item{sigma}{residual standard error.}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{df.residuals}{the residual degrees of freedom.}
  \item{formula}{\code{formula1} argument of \code{snp.quant}.}
  \item{rsquared}{adjusted r-squared values for the fitted model.}
  \item{AIC}{Akaike information criterion.}
  
}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A. McCaskie}
\note{

}
\seealso{\code{\link{epi.quant}}}
\examples{

data(pheno.dat)
mymodel <- epi.quant(formula=HDL~AGE+SBP, pheno=pheno.dat)
summary(mymodel)

# example with a subsetting variable
mymodel <- epi.quant(formula=HDL~AGE+SBP, pheno=pheno.dat, 
	sub=expression(SEX==1))
summary(mymodel)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

