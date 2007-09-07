\name{snp.cc.match}
\alias{snp.cc.match}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Single SNP analysis for matched case-control data}
\description{
\code{snp.cc.match} is used to fit conditional logistic regression models to single SNP genotype and phenotype, matched case-control data.  
}
\usage{
snp.cc.match(formula1, formula2, geno, pheno, sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula1}{a symbolic description of the full model to be fit, including SNP parameters. The response must be binary indicator of case-control status, and the formula must contain a variable indicating strata, or the matching sequence.}
  \item{formula2}{a symbolic description of the nested model excluding SNP parameters, to be compared to \code{formula1} in a likelihood ratio test. The response must be binary indicator of case-control status, and the formula must contain a variable indicating strata, or the matching sequence.}
  \item{geno}{a dataframe containing genotype data.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{sub}{an expression representing a subset of the data on which to perform the models.}
}
\details{
\code{formula1} should be in the form: \preformatted{response ~ predictor(s) + strata(strata_variable) + SNP(s)} and \code{formula2} should be in the form: \preformatted{response ~ predictor(s) + strata(strata_variable).} If case-control data is not matched, the \code{snp.bin} function should be used.

}
\value{
\code{snp.Clogit} returns an object of class \code{snpClogit}.

     The \code{summary} function can be used to obtain and print a
     summary of the results.  

     An object of class \code{snpClogit} is a list containing the
     following components:
     
  \item{results}{a table containing the hazard ratios, confidence intervals and p-values of the parameter estimates.}
  \item{formula}{\code{formula1} passed to \code{snp.cc.match}.}
  \item{Wald}{The Wald test for overall significance of the fitted model.}
  \item{logLik}{the log-likelihood for the model fit using \code{formula1}.}
  \item{fit.clogit}{an object of class \code{clogit} fit using \code{formula1}. See \code{clogit} for details.}
  \item{rsquared}{r-squared values for models fit using \code{formula1} and \code{formula2}.}
}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A. McCaskie}
\note{
}
\seealso{\code{\link{snp.bin}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by snp.cc.match
geno.dat <- SNP2Geno(SNP.dat, baseline=c("MM", "11", "GG", "CC"))

data(pheno.dat)
mymodel <- snp.cc.match(formula1=DISEASE~SBP+DBP+SNP_1_add+strata(STRAT), 
	formula2=DISEASE~SBP+DBP+strata(STRAT), pheno=pheno.dat, 
	geno=geno.dat)
summary(mymodel)

# example with subsetting variable
mymodel <- snp.cc.match(formula1=DISEASE~SBP+DBP+SNP_1_add+strata(STRAT), 
	formula2=DISEASE~SBP+DBP+strata(STRAT), pheno=pheno.dat, 
	geno=geno.dat, sub=expression(SEX==1))
summary(mymodel)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

