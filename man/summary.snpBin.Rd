\name{summary.snpBin}
\alias{summary.snpBin}
\alias{print.summary.snpBin}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing Single SNP analysis models}
\description{
  Summary method for objects of class \code{snpBin}
}
\usage{
\method{summary}{snpBin}(object, ...)
\method{print}{summary.snpBin}(x, digits = max(3, getOption("digits") - 3), 
	signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{snpBin}, the result of a call to \code{snp.bin}.}
  \item{x}{an object of class \code{summary.snpBin}, the result of a call to \preformatted{summary.snpBin.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
\code{summary.snpBin} returns an object of class \code{summary.snpBin}, a list
     with components

  \item{call}{the formula call.}
  \item{terms}{terms attribute of the formula called in \code{snp.bin}.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{df}{degrees of freedom parameter used in printing the model summary.}
  \item{residuals}{as per \code{glm}. The working residuals, that is the residuals in the final iteration of the IWLS fit.  Since cases with zero weights are omitted, their working residuals are \code{NA}.}
  \item{coefficients}{summarized results from fitted model, including odds ratios, confidence intervals and p-values.}
  \item{formula}{formula1 used in \code{snp.bin}.}
  \item{LRT}{a likelihood ratio test, testing for significant improvement of the model when SNP parameters are included.}
  \item{weights}{as per \code{glm}. The working weights, that is the weights in the final iteration of the IWLS fit.}
  \item{AIC}{Akaike Information Criterion for the generalized linear model fit including SNPs.}

}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A. McCaskie}

\seealso{\code{\link{snp.bin}}, \code{\link{snp.quant}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by snp.bin
geno.dat <- SNP2Geno(SNP.dat, baseline=c("MM", "11", "GG", "CC"))

data(pheno.dat)
mymodel <- snp.bin(formula1=PLAQUE~AGE+SEX+SNP_1_add, 
	formula2=PLAQUE~AGE+SEX, geno=geno.dat, pheno=pheno.dat)
summary(mymodel)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{methods}
