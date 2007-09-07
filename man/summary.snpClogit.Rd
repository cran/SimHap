\name{summary.snpClogit}
\alias{summary.snpClogit}
\alias{print.summary.snpClogit}


\title{Summarizing Single SNP analysis models for matched case-control data}
\description{
  Summary method for objects of class \code{snpClogit}
}
\usage{
\method{summary}{snpClogit}(object, ...)
\method{print}{summary.snpClogit}(x, digits = max(3, getOption("digits") - 3), 
	signif.stars = getOption("show.signif.stars"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{snpClogit}}
  \item{x}{an object of class \code{summary.snpClogit}, the result of a call to \preformatted{summary.snpClogit.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.snpClogit} returns an object of class \code{summary.snpClogit}, a list with components
  \item{terms}{terms attribute of formula1 called in \code{snp.cc.match}.}
  \item{coefficients}{summarized results from fitted model, including odds ratios and p-values.}
  \item{formula}{formula1 used in \code{snp.cc.match}.}
  \item{LRT}{likelihood ratio test comparing the model with SNP variables compared to the model without SNPs.}
  \item{Wald}{Wald statistic for the fitted model.}
  \item{rsquared}{adjusted r-squared values for the fitted model.}
  \item{residuals}{the residuals, that is response minus fitted values.}

}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A. McCaskie}
\seealso{\code{\link{snp.cc.match}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by snp.cc.match
geno.dat <- SNP2Geno(SNP.dat, baseline=c("MM", "11", "GG", "CC"))

data(pheno.dat)
mymodel <- snp.cc.match(formula1=DISEASE~SBP+DBP+SNP_1_add+strata(STRAT), 
	formula2=DISEASE~SBP+DBP+strata(STRAT), pheno=pheno.dat, 
	geno=geno.dat)
summary(mymodel)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

