\name{summary.snpSurv}
\alias{summary.snpSurv}
\alias{print.summary.snpSurv}

\title{Summarizing Single SNP analysis models for survival data}
\description{
  Summary method for objects of class \code{snpSurv}
}
\usage{
\method{summary}{snpSurv}(object, ...)
\method{print}{summary.snpSurv}(x, digits = max(3, getOption("digits") - 3), 
	signif.stars = getOption("show.signif.stars"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{snpSurv}}
  \item{x}{an object of class \code{summary.snpSurv}, the result of a call to \preformatted{summary.snpSurv.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.snpSurv} returns an object of class \code{summary.snpSurv}, a list with components
  \item{terms}{terms attribute of formula1 called in \code{snp.surv}.}
  \item{coefficients}{summarized results from fitted model, including hazard ratios and p-values.}
  \item{formula}{formula1 used in \code{snp.surv}.}
  \item{LRT}{likelihood ratio test comparing the model with SNP variables compared to the model without SNPs.}
  \item{Wald}{Wald statistic for the fitted model.}
  \item{rsquared}{adjusted r-squared values for the fitted model.}
  \item{residuals}{the residuals, that is response minus fitted values.}

}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A. McCaskie}
\seealso{\code{\link{snp.surv}}}
\examples{

data(SNPsurv.dat)

# transforms SNPlong.dat to an object containing 3 columns
# per SNP - additive, dominant and recessive, where genotypes
# defined in 'baseline' serve as the baseline genotypes
survGeno.dat <- SNP2Geno(SNPsurv.dat, baseline=c("V2V2", "GG", "CC"))

data(survPheno.dat)
mymodel <- snp.surv(formula1=Surv(time, status)~age+SNP_1_add, 
	formula2=Surv(time, status)~age, pheno=survPheno.dat, 
	geno=survGeno.dat)
summary(mymodel)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

