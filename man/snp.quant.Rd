\name{snp.quant}
\alias{snp.quant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Single SNP analysis for quantitative outcomes}
\description{
\code{snp.quant} is used to fit linear regression models to single SNP genotype and phenotype data for continuous Normal outcomes. 
}
\usage{
snp.quant(formula1, formula2, geno, pheno, sub = NULL, 
	predict_variable = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula1}{a symbolic description of the full model to be fit, including SNP parameters. The details of model specification are given below.}
  \item{formula2}{a symbolic description of the nested model excluding SNP parameters, to be compared to \code{formula1} in a likelihood ratio test.}
  \item{geno}{a dataframe containing genotype data.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{sub}{an expression representing a subset of the data on which to perform the models.}
  \item{predict_variable}{a factored SNP variable, at each level of which the response variable will be predicted based on the model fit using formula1.}
}
\details{
\code{formula1} should be in the form: \code{response ~ predictor(s) + SNP(s)} and \code{formula2} should be in the form: \code{response ~ predictor(s)}. A formula has an implied intercept term. See documentation for \code{formula} function for more details of allowed formulae.

}
\value{
 \code{snp.quant} returns an object of 'class' \code{snpQuant}.

     The \code{summary} function can be used to obtain and print a
     summary of the results.  

     An object of class \code{snpQuant} is a list containing the
     following components:
     
  \item{results}{a table containing the coefficients, standard errors and p-values of the parameter estimates.}
  \item{formula1}{\code{formula1} passed to \code{snp.quant}.}
  \item{formula2}{\code{formula2} passed to \code{snp.quant}.}
  \item{LRT}{a likelihood ratio test, testing for significant improvement of the model when genotypic parameters are included.}
  \item{ANOD}{analysis of deviance table for the model fit using \code{formula1}.}
  \item{logLik}{the log-likelihood for the linear model fit using \code{formula1}.}
  \item{fit.lm}{a \code{lm} object fit using \code{formula1}.}
  \item{fitsub.lm}{a \code{lm} object fit using \code{formula2}.}
  \item{rsquared}{r-squared values for models fit using \code{formula1} and \code{formula2}.}
  \item{predicted.values}{estimated marginal means of the outcome variable broken down by SNP levels, evaluated at mean values of the model predictors.}
  \item{AIC}{Akaike Information Criterion for the linear model fit using \code{formula1}.}

}
\references{
Chambers, J.M. (1992) Linear models. Chapter 4 of \emph{Statistical Models in S}, eds Chambers, J.M., Hastie, T.J., Wadsworth & Brooks/Cole.

McCaskie, P.A., Carter, K.W, Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

Wilkinson, G.N., Rogers, C.E. (1973) Symbolic descriptions of factorial models for analysis of variance. \emph{Applied Statistics}, 22, 392-9.

}
\author{Pamela A. McCaskie}

\seealso{\code{\link{snp.bin}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by snp.quant
geno.dat <- SNP2Geno(SNP.dat, baseline=c("MM", "11", "GG", "CC"))

data(pheno.dat)
mymodel <- snp.quant(formula1=HDL~AGE+SBP+factor(SNP_1_add),
	formula2=HDL~AGE+SBP, geno=geno.dat, pheno=pheno.dat)
summary(mymodel)

# example using a variable for which to predict marginal means
mymodel <- snp.quant(formula1=HDL~AGE+SBP+factor(SNP_1_add),
	formula2=HDL~AGE+SBP, geno=geno.dat, pheno=pheno.dat, 
	predict_variable="SNP_1_add")
summary(mymodel)

# example with a subsetting variable, looking at males only
mymodel <- snp.quant(formula1=HDL~AGE+SBP+factor(SNP_1_add),
	formula2=HDL~AGE+SBP, geno=geno.dat, pheno=pheno.dat, 
	sub=expression(SEX==1))
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

