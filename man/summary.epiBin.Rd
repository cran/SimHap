\name{summary.epiBin}
\alias{summary.epiBin}
\alias{print.summary.epiBin}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing epidemiological analysis models for binary outcomes}
\description{
  Summary method for objects of class \code{epiBin}
}
\usage{
\method{summary}{epiBin}(object, ...)
\method{print}{summary.epiBin}(x, digits = max(3, getOption("digits") - 3),
	signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{epiBin}, a result of a call to \code{epi.bin}.}
  \item{x}{an object of class \code{summary.epiBin}, the result of a call to \code{summary.epiBin}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.epiBin} returns an object of class \code{summary.epiBin}, a list
     with components

  \item{call}{the formula call.}
  \item{terms}{terms attribute of the formula called in \code{epi.bin}.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{df}{degrees of freedom parameter used in printing the model summary.}
  \item{coefficients}{summarized results from fitted model, including odds ratios, confidence intervals and p-values.}
  \item{residuals}{as per \code{glm}. The working residuals, that is the residuals in the final iteration of the IWLS fit.  Since cases with zero weights are omitted, their working residuals are \code{NA}.}
  \item{formula}{formula1 used in \code{epi.bin}.}
  \item{weights}{as per \code{glm}. The working weights, that is the weights in the final iteration of the IWLS fit.}
  \item{AIC}{Akaike Information Criterion for the generalized linear model fit including SNPs.}

}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

}
\author{Pamela A. McCaskie}

\seealso{\code{\link{epi.bin}}}
\examples{

data(pheno.dat)
mymodel <- epi.bin(formula=PLAQUE~AGE+SEX, pheno=pheno.dat)
summary(mymodel)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

