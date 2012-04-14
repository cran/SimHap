\name{summary.snpQuant}
\alias{summary.snpQuant}
\alias{print.summary.snpQuant}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing single SNP analysis models for quantitative outcomes}
\description{
Summary method for objects of class \code{snpQuant}
}
\usage{
\method{summary}{snpQuant}(object, ...)
\method{print}{summary.snpQuant}(x, digits = max(3, getOption("digits") - 3), 
	signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{object}{an object of class \code{snpQuant}, a result of a call to \code{snp.quant}.}
  \item{x}{an object of class \code{summary.snpQuant}, the result of a call to \preformatted{summary.snpQuant.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  \code{summary.snpQuant} returns an object of class \code{summary.snpQuant}, a list with components
  
  \item{call}{the formula call.}
  \item{terms}{terms attribute of the formula called in \code{snp.quant}.}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{df.residuals}{the residual degrees of freedom.}
  \item{na.action}{method used for missing data.}
  \item{df}{residual degrees of freedom.}
  \item{sigma}{residual standard error.}
  \item{coefficients}{summarized results from fitted model, including coefficients, standard errors and p-values.}
  \item{formula}{formula1 used in \code{snp.quant}.}
  \item{LRT}{likelihood ratio test comparing the model with SNP variables compared to the model without SNPs.}
  \item{AIC}{Akaike information criterion for the full fitted model.}
  \item{rsquared}{adjusted r-squared values for the fitted model.}
  \item{predicted.values}{estimated marginal means for a chosen model term.}
}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A. McCaskie}

\seealso{\code{\link{snp.quant}}, \code{\link{snp.bin}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by snp.quant
geno.dat <- SNP2Geno(SNP.dat, baseline=c("MM", "11", "GG", "CC"))

data(pheno.dat)
mymodel <- snp.quant(formula1=HDL~AGE+SBP+factor(SNP_1_add),
	formula2=HDL~AGE+SBP, geno=geno.dat, pheno=pheno.dat)
summary(mymodel)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
