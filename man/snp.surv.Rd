\name{snp.surv}
\alias{snp.surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Single SNP analysis for survival data}
\description{
\code{snp.surv} is used to fit cox proportional hazards models to single SNP genotype and phenotype survival data.  
}
\usage{
snp.surv(formula1, formula2, geno, pheno, sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula1}{a symbolic description of the full model to be fit, including SNP parameters. The response must be a survival object as returned by the \code{Surv} function.}
  \item{formula2}{a symbolic description of the nested model excluding SNP parameters, to be compared to \code{formula1} in a likelihood ratio test. The response must be a survival object as returned by the \code{Surv} function.}
  \item{geno}{a dataframe containing genotype data.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{sub}{an expression representing a subset of the data on which to perform the models.}
}
\details{
\code{formula1} should be in the form: \code{response ~ predictor(s) + SNP(s)} and \code{formula2} should be in the form: \code{response ~ predictor(s)}. A formula has an implied intercept term. See documentation for the \code{formula} function for more details of allowed formulae.

}
\value{
\code{snp.surv} returns an object of class \code{snpSurv}.

     The \code{summary} function can be used to obtain and print a
     summary of the results.  

     An object of class \code{snpSurv} is a list containing the
     following components:
     
  \item{results}{a table containing the hazard ratios, confidence intervals and p-values of the parameter estimates.}
  \item{formula}{\code{formula1} passed to \code{snp.surv}.}
  \item{Wald}{The Wald test for overall significance of the fitted model including SNP parameters.}
  \item{logLik}{the log-likelihood for the model fit using \code{formula1}.}
  \item{fit.coxph}{an object of class \code{coxph} fit using \code{formula1}. See \code{coxph.object} for details.}
  \item{rsquared}{r-squared values for models fit using \code{formula1} and \code{formula2}.}
}
\references{
Andersen, P., Gill, R. (1982) Cox's regression model for counting processes, a large sample study, \emph{Annals of Statistics}, 10:1100-1120.

McCaskie, P.A., Carter, K.W, Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

Therneau, T., Grambsch, P., Fleming, T. Martingale based residuals for survival models,  \emph{Biometrika}, 77(1):147-160.
}
\author{Pamela A. McCaskie}

\seealso{\code{\link{haplo.surv}}}
\examples{

data(SNPsurv.dat)

# transforms SNPlong.dat to an object containing 3 columns
# per SNP - additive, dominant and recessive, where genotypes
# defined in 'baseline' serve as the baseline genotypes
survGeno.dat <- SNP2Geno(SNPsurv.dat, baseline=c("V2V2", "GG", "CC"))

data(survPheno.dat)
mymodel <- snp.surv(formula1=Surv(time, status)~age+SNP_1_add, 
	formula2=Surv(time, status)~age, geno=survGeno.dat, 
	pheno=survPheno.dat)
summary(mymodel)

# example with subsetting variable
mymodel <- snp.surv(formula1=Surv(time, status)~age+SNP_1_add, 
	formula2=Surv(time, status)~age, pheno=survPheno.dat, 
	geno=survGeno.dat, sub=expression(sex==1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

