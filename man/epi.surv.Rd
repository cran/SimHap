\name{epi.surv}
\alias{epi.surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Epidemiological analysis for survival data}
\description{
\code{epi.surv} is used to fit Cox proportional hazards models to epidemiological survival data.
}
\usage{
epi.surv(formula, pheno, sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The response must be a survival object as returned by the \code{Surv} function.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{sub}{an expression representing a subset of the data on which to perform the models.}
}
\details{
\code{formula} should be in the form of \code{response ~ predictor(s)}. A formula has an implied intercept term. See documentation for the \code{formula} function for more details of allowed formulae.
}
\value{
  \code{epi.surv} returns an object of class \code{epiSurv} containing the following items
  
  \item{results}{a table containing the hazard ratios, confidence intervals and p-values of the parameter estimates.}
  \item{formula}{\code{formula} passed to \code{epi.surv}.}
  \item{Wald}{The Wald test for overall significance of the fitted model including SNP parameters.}
  \item{logLik}{the log-likelihood for the linear model fit using \code{formula}.}
  \item{fit.coxph}{an object of class \code{coxph} fit using \code{formula1}. See \code{coxph.object} for details.}
  \item{rsquared}{r-squared values for the fitted model.}

}
\references{
Andersen, P., Gill, R. (1982) Cox's regression model for counting processes, a large sample study, \emph{Annals of Statistics}, 10:1100-1120.

McCaskie, P.A., Carter, K.W, Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

Therneau, T., Grambsch, P., Fleming, T. Martingale based residuals for survival models,  \emph{Biometrika}, 77(1):147-160.
}
\author{Pamela A McCaskie}

\seealso{\code{\link{snp.surv}}, \code{\link{haplo.surv}}}
\examples{

data(survPheno.dat)
mymodel <- epi.surv(formula=Surv(time, status)~age, pheno=survPheno.dat)
summary(mymodel)

# example with subsetting variable
mymodel <- epi.surv(formula=Surv(time, status)~age, pheno=survPheno.dat, 
	sub=expression(sex==1))
summary(mymodel)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

