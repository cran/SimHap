\name{snp.bin}
\alias{snp.bin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Single SNP analysis for binary outcomes}
\description{
  \code{snp.bin} is used to fit generalized linear regression models to single SNP genotype and phenotype data for a binary outcome.
}
\usage{
snp.bin(formula1, formula2, geno, pheno, sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula1}{a symbolic description of the full model to be fit, including SNP parameters. The details of model specification are given below.}
  \item{formula2}{a symbolic description of the nested model excluding SNP parameters, to be compared to \code{formula1} in a likelihood ratio test.}
  \item{geno}{a dataframe containing genotype data.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{sub}{an expression representing a subset of the data on which to perform the models.}
}
\details{
\code{formula1} should be in the form of \code{outcome ~ predictor(s) + SNP(s)} and \code{formula2} should be in the form \code{outcome ~ predictor(s)}. A formula has an implied intercept term. See documentation for \code{formula} function for more details of allowed formulae.
}
\value{

\code{snp.bin} returns an object of 'class' \code{snpBin}.

     The \code{summary} function can be used to obtain and print a
     summary of the results.  

     An object of class \code{snpBin} is a list containing the
     following components:
     
  \item{results}{a table containing the odds ratios, confidence intervals and p-values of the parameter estimates.}
  \item{formula1}{formula1 passed to snp.bin.}
  \item{formula2}{formula2 passed to snp.bin.}
  \item{LRT}{a likelihood ratio test, testing for significant improvement of the model when haplotypic parameters are included.}
  \item{ANOD}{analysis of deviance table for the model fit using formula1.}
  \item{logLik}{the log-likelihood for the linear model fit using formula1.}
  \item{fit.glm}{a \code{glm} object fit using formula1.}
  \item{fitsub.glm}{a \code{glm} object fit using formula2.}
  \item{AIC}{Akaike Information Criterion for the linear model fit using formula1.}
  
}
\references{
Dobson, A.J. (1990) \emph{An Introduction to Generalized Linear Models}. London: Chapman and Hall.

Hastie, T.J., Pregibon, D. (1992) Generalized linear models. Chapter 6 of \emph{Statistical Models in S}, eds Chambers, J.M., Hastie, T.J., Wadsworth & Brooks/Cole.

McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

McCullagh, P., Nelder, J.A. (1989) \emph{Generalized Linear Models}. London: Chapman and Hall.

Venables, W.N., Ripley, D.B. (2002) \emph{Modern Applied Statistics with S}. New York: Springer.
}
\author{Pamela A. McCaskie}

\seealso{\code{\link{snp.quant}}, \code{\link{haplo.quant}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by snp.bin
geno.dat <- SNP2Geno(SNP.dat, baseline=c("MM", "11", "GG", "CC"))
data(pheno.dat)

mymodel <- snp.bin(formula1=PLAQUE~AGE+SEX+SNP_1_add, 
	formula2=PLAQUE~AGE+SEX, geno=geno.dat, pheno=pheno.dat)
summary(mymodel)

# example with a subsetting variable, looking at 
# people over 50 years of age only
mymodel <- snp.bin(formula1=PLAQUE~AGE+SEX+SNP_1_add, 
	formula2=PLAQUE~AGE+SEX, geno=geno.dat, pheno=pheno.dat, 
	sub=expression(AGE>50))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
