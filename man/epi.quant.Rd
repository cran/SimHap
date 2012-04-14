\name{epi.quant}
\alias{epi.quant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Epidemiological analysis for quantitative outcomes}
\description{
\code{epi.quant} is used to fit linear regression models to single SNP genotype and phenotype data for a continuous Normal outcome. 
}
\usage{
epi.quant(formula, pheno, sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the full model to be fit. The details of model specification are given below.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{sub}{an expression representing a subset of the data on which to perform the models.}
}
\details{
\code{formula} should be in the form of \code{response ~ predictor(s)}. A formula has an implied intercept term. See documentation for \code{formula} function for more details of allowed formulae.

}
\value{
  \code{epi.quant} returns an object of class \code{epiQuant} containing the following items
  
  \item{formula}{formula passed to epi.quant.}
  \item{results}{a table containing the coefficients, standard errors and p-values of the parameter estimates.}
  \item{fit.lm}{a \code{lm} object fit using formula.}
   \item{ANOD}{analysis of deviance table for the model fit using formula.}
  \item{logLik}{the log-likelihood for the linear model fit using formula.}
  \item{AIC}{Akaike Information Criterion for the linear model fit using formula.}

}
\references{
Chambers, J.M. (1992) Linear models. Chapter 4 of \emph{Statistical Models in S}, eds Chambers, J.M., Hastie, T.J., Wadsworth & Brooks/Cole.

McCaskie, P.A., Carter, K.W, Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

Wilkinson, G.N., Rogers, C.E. (1973) Symbolic descriptions of factorial models for analysis of variance. \emph{Applied Statistics}, 22, 392-9.

}
\author{Pamela A. McCaskie}

\seealso{\code{\link{snp.quant}}, \code{\link{haplo.quant}}, \code{\link{epi.bin}}}
\examples{

data(pheno.dat)
mymodel <- epi.quant(formula=LDL~AGE+SBP, pheno=pheno.dat)
summary(mymodel)

# example with a subsetting variable, looking at males only
mymodel <- epi.quant(formula=LDL~AGE+SBP, pheno=pheno.dat, 
	sub=expression(SEX==1))
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
