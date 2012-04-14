\name{summary.epiClogit}
\alias{summary.epiClogit}
\alias{print.summary.epiClogit}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing epidemiological analysis models for matched case-control data}
\description{
  Summary method for objects of class \code{epiClogit}
}
\usage{
\method{summary}{epiClogit}(object, ...)
\method{print}{summary.epiClogit}(x, digits = max(3, getOption("digits") - 3), 
	signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{object}{an object of class \code{epiClogit}, a result of a call to \code{epi.cc.match}.}
  \item{x}{an object of class \code{summary.epiClogit}, the result of a call to \preformatted{summary.epiClogit.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.epiClogit} returns an object of class summary.snpClogit, a list with components
  
  \item{terms}{terms attribute of formula called in \code{epi.cc.match}.}
  \item{coefficients}{summarized results from fitted model, including odds ratios and p-values.}
  \item{formula}{formula used in \code{epi.Clogit}.}
  \item{Wald}{Wald statistic for the fitted model.}
  \item{rsquared}{adjusted r-squared values for the fitted model.}
  \item{residuals}{the residuals, that is response minus fitted values.}
}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

}
\author{Pamela A McCaskie}

\seealso{\code{\link{epi.cc.match}}}
\examples{

data(pheno.dat)
mymodel <- epi.cc.match(formula=DISEASE~SBP+DBP+strata(STRAT), 
	pheno=pheno.dat)
summary(mymodel)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

