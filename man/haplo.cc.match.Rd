\name{haplo.cc.match}
\alias{haplo.cc.match}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Haplotype analysis for matched case-control data}
\description{
\code{haplo.surv} performs a series of conditional logistic regression models to matched case-control data with haplotypes using a simulation-based approach to account for uncertainty in haplotype assignment when phase is unknown. 
}
\usage{
haplo.cc.match(formula1, formula2, pheno, haplo, sim, effect = "add", 
	sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula1}{a symbolic description of the full model to be fit, including haplotype parameters. The response must be binary indicator of case-control status, and the formula must contain a variable indicating strata, or the matching sequence.}
  \item{formula2}{a symbolic description of the nested model excluding haplotype parameters, to be compared to \code{formula1} in a likelihood ratio test. The response must be binary indicator of case-control status, and the formula must contain a variable indicating strata, or the matching sequence.}
  \item{pheno}{a dataframe containing phenotype data.}
  \item{haplo}{a haplotype object made by \code{make.haplo.rare}.}
  \item{sim}{the number of simulations from which to evaluate the results.}
  \item{effect}{the genetic effect type: \code{"add"} for additive, \code{"dom"} for dominant and \code{"rec"} for recessive. Defaults to additive. See note.}
  \item{sub}{optional. An expression using a binary operator, representing a subset of individuals on which to perform analysis. e.g. \code{sub=expression(sex==1)}.}
}
\details{
  \code{formula1} should be in the form: \preformatted{response ~ predictor(s) + strata(strata_variable) + haplotype(s)} and \code{formula2} should be in the form: \preformatted{response ~ predictor(s) + strata(strata_variable)}. If case-control data is not matched, the \code{haplo.bin} function should be used. 
}
\value{
\code{haplo.cc.match} returns an object of 'class' \code{hapClogit}.

     The \code{summary} function can be used to obtain and print a
     summary of the results.  

     An object of class \code{hapClogit} is a list containing the
     following components:
     
  \item{formula1}{\code{formula1} passed to \code{haplo.cc.match}.}
  \item{formula2}{\code{formula2} passed to \code{haplo.cc.match}.}
  \item{results}{a table containing the odds ratios, confidence intervals and p-values of the parameter estimates, averaged over the n=\code{sim} models performed.}
  \item{empiricalResults}{a list containing the odds ratios, confidence intervals and p-values calculated at each simulation}
  \item{logLik}{the average log-likelihood for the n=\code{sim} linear models fit using \code{formula1}.}
  \item{LRT}{a likelihood ratio test, testing for significant improvement of the model when haplotypic parameters are included}
  \item{ANOVA}{analysis of variance, comparing the two models fit with and without haplotypic parameters.}
  \item{Wald}{The Wald test for overall significance of the fitted model including haplotypes.}
  \item{rsquared}{r-squared values for models fit using \code{formula1} and \code{formula2}.}
  \item{effect}{the haplotypic effect modelled, `ADDITIVE', `DOMINANT' or `RECESSIVE'.}
}
\references{
Little, R.J.A., Rubin, D.B. (2002) \emph{Statistical Analysis with Missing Data}. John Wiley and Sons, New Jersey. 

McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

Rubin, D.B. (1996) Multiple imputation after 18+ years (with discussion). \emph{Journal of the American Statistical Society}, 91:473-489.
}
\author{Pamela A. McCaskie}
\note{
To model a codominant haplotypic effect, define the desired haplotype as a factor in the \code{formula1} argument. e.g. \code{factor(h.AAA)}, and use the default option for \code{effect}.
}
\seealso{\code{\link{snp.cc.match}}, \code{\link{haplo.bin}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by infer.haplos
haplo.dat <- SNP2Haplo(SNP.dat)

data(pheno.dat)

# generate haplotype frequencies and haplotype design matrix
myinfer<-infer.haplos(haplo.dat) 

# prints haplotype frequencies among cases
myinfer$hap.freq.cases

# prints haplotype frequencies among controls
myinfer$hap.freq.controls

# generate haplo object where haplotypes with a frequency 
# below min.freq are grouped as a category called "rare"
myhaplo<-make.haplo.rare(myinfer,min.freq=0.05)

mymodel <- haplo.cc.match(formula1=DISEASE~SBP+DBP+h.N1AA+strata(STRAT), 
	formula2=DISEASE~SBP+DBP+strata(STRAT), haplo=myhaplo, 
	pheno=pheno.dat, sim=10)
summary(mymodel)

# example using a subsetting variable - looking at males only
mymodel <- haplo.cc.match(formula1=DISEASE~SBP+DBP+h.N1AA+strata(STRAT), 
	formula2=DISEASE~SBP+DBP+strata(STRAT), haplo=myhaplo, 
	pheno=pheno.dat, sim=10, sub=expression(SEX==1))
summary(mymodel)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

