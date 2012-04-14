\name{epi.bin}
\alias{epi.bin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Epidemiological analysis for binary outcomes}
\description{
\code{epi.bin} is used to fit generalized linear regression models to epidemiological phenotype data for a binary outcome, assuming a binomial error distribution and logit link function.
}
\usage{
epi.bin(formula, pheno, sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details of model specification are given below.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{sub}{an expression representing a subset of the data on which to perform the model.}
}
\details{
\code{formula} should be in the form of \code{outcome ~ predictor(s)}. A formula has an implied intercept term. See documentation for \code{formula} function for more details of allowed formulae.
}
\value{
\code{epi.bin} returns an object of class \code{epiBin} containing the following items
  
  \item{formula}{formula passed to \code{epi.bin}.}
  \item{results}{a table containing the odds ratios, confidence intervals and p-values of the parameter estimates.}
  \item{fit.glm}{a \code{glm} object fit using \code{formula}.}
  \item{ANOD}{analysis of deviance table for the model fit using \code{formula}.}
  \item{logLik}{the log-likelihood for the linear model fit using \code{formula}.}
  \item{AIC}{Akaike Information Criterion for the linear model fit using \code{formula}.}
 
}
\references{
Dobson, A.J. (1990) \emph{An Introduction to Generalized Linear Models}. London: Chapman and Hall.

Hastie, T.J., Pregibon, D. (1992) Generalized linear models. Chapter 6 of \emph{Statistical Models in S}, eds Chambers, J.M., Hastie, T.J., Wadsworth & Brooks/Cole.

McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

McCullagh, P., Nelder, J.A. (1989) \emph{Generalized Linear Models}. London: Chapman and Hall.

Venables, W.N., Ripley, D.B. (2002) \emph{Modern Applied Statistics with S}. New York: Springer.
}
\author{Pamela A. McCaskie}

\seealso{\code{\link{snp.bin}}, \code{\link{haplo.bin}}, \code{\link{epi.quant}}}
\examples{

data(pheno.dat)
mymodel <- epi.bin(formula=PLAQUE~AGE+SBP, pheno=pheno.dat)
summary(mymodel)

# example with a subsetting variable, looking at males only
mymodel <- epi.bin(formula=PLAQUE~AGE+SBP, pheno=pheno.dat, 
	sub=expression(SEX==1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

