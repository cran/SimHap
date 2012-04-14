\name{epi.cc.match}
\alias{epi.cc.match}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Epidemiological analysis for matched case-control data}
\description{
\code{epi.cc.match} is used to fit conditional logistic regression models to matched case-control data.
}
\usage{
epi.cc.match(formula, pheno, sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The response must be binary indicator of case-control status, and the formula must contain a variable indicating strata, or the matching sequence.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{sub}{an expression representing a subset of the data on which to perform the models.}
}
\details{
\code{formula} should be in the form: \preformatted{response ~ predictor(s) + strata(strata_variable).} 
}
\value{
  \code{epi.cc.match} returns an object of class \code{epiClogit} containing the following items
  
  \item{results}{a table containing the odds ratios, confidence intervals and p-values of the parameter estimates.}
  \item{formula}{\code{formula} passed to \code{epi.cc.match}.}
  \item{Wald}{The Wald test for overall significance of the fitted model including SNP parameters.}
  \item{logLik}{the log-likelihood for the linear model fit using \code{formula}.}
  \item{fit.coxph}{an object of class \code{clogit} fit using \code{formula1}. See \code{clogit} for details.}
  \item{rsquared}{r-squared values for the fitted model.}

}
\references{

McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A McCaskie}

\seealso{\code{\link{epi.bin}}, \code{\link{clogit}}}
\examples{

data(pheno.dat)
mymodel <- epi.cc.match(formula=DISEASE~SBP+DBP+strata(STRAT), 
	pheno=pheno.dat)
summary(mymodel)

# example with subsetting variable
mymodel <- epi.cc.match(formula=DISEASE~SBP+DBP+strata(STRAT), 
	pheno=pheno.dat, sub=expression(SEX==1))
summary(mymodel)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

