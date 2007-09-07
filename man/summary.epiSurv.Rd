\name{summary.epiSurv}
\alias{summary.epiSurv}
\alias{print.summary.epiSurv}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing epidemiological analysis models for survival data}
\description{
  Summary method for objects of class \code{epiSurv}
}
\usage{
\method{summary}{epiSurv}(object, ...)
\method{print}{summary.epiSurv}(x, digits = max(3, getOption("digits") - 3), 
	signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{object}{an object of class \code{epiSurv}, a result of a call to \code{epi.surv}.}
  \item{x}{an object of class \code{summary.epiSurv}, the result of a call to \preformatted{summary.epiSurv.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.epiSurv} returns an object of class summary.epiSurv, a list with components
  
  \item{terms}{terms attribute of formula called in \code{epi.surv}.}
  \item{coefficients}{summarized results from fitted model, including hazard ratios and p-values.}
  \item{formula}{formula used in \code{epi.surv}.}
  \item{Wald}{Wald statistic for the fitted model.}
  \item{rsquared}{adjusted r-squared values for the fitted model.}
  \item{residuals}{the residuals, that is response minus fitted values.}
}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A McCaskie}

\seealso{\code{\link{epi.surv}}}
\examples{

data(survPheno.dat)
mymodel <- epi.surv(formula=Surv(time, status)~age, pheno=survPheno.dat)
summary(mymodel)

# example with subsetting variable
mymodel <- epi.surv(formula=Surv(time, status)~age, pheno=survPheno.dat, 
	sub=expression(sex==1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

