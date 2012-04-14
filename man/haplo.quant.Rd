\name{haplo.quant}
\alias{haplo.quant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Haplotype analysis for a Normally distributed quantitative trait}
\description{
  \code{haplo.quant} performs a series of linear models using a simulation-based approach to account for uncertainty in haplotype assignment when phase is unknown. 
}
\usage{
haplo.quant(formula1, formula2, pheno, haplo, sim, 
	effect = "add", sub = NULL, predict_variable = NULL, adjust = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula1}{a symbolic description of the full model including haplotype parameters to be fit. The details of model specification are given below.}
  \item{formula2}{a symbolic description of the nested model excluding haplotype parameters, to be compared to \code{formula1} in a likelihood ratio test.}
  \item{pheno}{a dataframe containing phenotype data. }
  \item{haplo}{a haplotype object made by \code{make.haplo.rare}. The
    subjects must in the same order as they are in the phenotype data.}
  \item{sim}{the number of simulations from which to evaluate the results.}
  \item{effect}{the haplotypic effect type: \code{"add"} for additive, \code{"dom"} for dominant and \code{"rec"} for recessive. Defaults to additive. See note.}
  \item{sub}{optional. An expression representing a subset of individuals on which to perform analysis. e.g. \code{sub=expression(sex==1)}.}
  \item{predict_variable}{an expression using a binary operator,
    representing a subset of the data on which to perform the models}
  \item{adjust}{a logical flag. If adjust=TRUE, the adjusted degrees of
    freedom is used. This is recommended when the computed degrees of
    freedom is larger than the complete data degrees of freedom. By
    default, adjust=FALSE.}

  
}
\details{
  \code{formula1} should be in the form of \code{response ~ predictor(s) + haplotype(s)} and \code{formula2} should be in the form \code{response ~ predictor(s)}. A formula has an implied intercept term. See 'formula' for more details of allowed formulae.
}
\value{
    \code{haplo.quant} returns an object of class \code{hapQuant}.

     The \code{summary} function can be used to obtain and print a
     summary of the results.  

     An object of class \code{hapQuant} is a list containing the
     following components:

  \item{formula1}{\code{formula1} passed to \code{haplo.quant}.}
  \item{formula1}{\code{formula2} passed to \code{haplo.quant}.}
  \item{results}{a table containing the coefficients, averaged over the
    \code{sim} models performed; standard errors, computed as the sum of
    the between-imputation and within-imputation variance; and
    p-values, based on a t-distribution with appropriately computed
  degrees of freedom, of the parameter estimates.}
  \item{empiricalResults}{a list containing the coefficients, standard errors and p-values calculated at each simulation.}
  \item{rsquared}{r-squared values for models fit using \code{formula1} and \code{formula2}.}
  \item{ANOD}{analysis of deviance table for the model fit using \code{formula1}, averaged over all simulations.}
  \item{logLik}{the average log-likelihood for the linear model fit using \code{formula1}.}
  \item{WALD}{a Wald test, testing for significant improvement of the
    model when haplotypic parameters are included.}
  \item{predicted}{estimated marginal means of the outcome variable broken down by haplotype levels, evaluated at mean values of the model predictors, averaged over all simulations.}
  \item{empiricalPredicted}{estimated marginal means calculated at each simulation.}
  \item{aic}{Akaike Information Criterion for the linear model fit using formula1, averaged over all simulations.}
  \item{aicPredicted}{Akaike Information Criteria calculated at each simulation.}
  \item{effect}{the haplotypic effect modelled, `ADDITIVE', `DOMINANT' or `RECESSIVE'.}

}
\references{
Chambers, J.M. (1992) Linear models. Chapter 4 of \emph{Statistical Models in S}, eds Chambers, J.M., Hastie, T.J., Wadsworth & Brooks/Cole.

Little, R.J.A., Rubin, D.B. (2002) \emph{Statistical Analysis with Missing Data}. John Wiley and Sons, New Jersey. 

McCaskie, P.A., Carter, K.W, Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

Rubin, D.B. (1996) Multiple imputation after 18+ years (with discussion). \emph{Journal of the American Statistical Society}, 91:473-489.

Wilkinson, G.N., Rogers, C.E. (1973) Symbolic descriptions of factorial
models for analysis of variance. \emph{Applied Statistics}, 22, 392-9.

Barnard, J., Rubin, D.B. (1999) Small-sample degrees of freedom with
multiple imputation. \emph{Biometrika}, 86, 948-955.

}
\author{Pamela A. McCaskie}
\note{To model a codominant haplotypic effect, define the desired haplotype as a factor in the \code{formula1} argument. e.g. \code{factor(h.AAA)}, and use the default option for \code{effect}

}
\seealso{\code{\link{haplo.bin}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by infer.haplos
haplo.dat <- SNP2Haplo(SNP.dat)

data(pheno.dat)

# generate haplotype frequencies and haplotype design matrix
myinfer<-infer.haplos(haplo.dat)

# print haplotype frequencies generated by infer.haplos
myinfer$hap.freq 

# generate haplo object where haplotypes with a frequency 
# below min.freq are grouped as a category called "rare"
myhaplo<-make.haplo.rare(myinfer,min.freq=0.05) 
mymodel <- haplo.quant(formula1=HDL~AGE+SBP+h.N1AA, 
	formula2=HDL~AGE+SBP, pheno=pheno.dat, haplo=myhaplo, sim=10)
summary(mymodel)

# example using a variable for which to predict marginal means
mymodel <- haplo.quant(formula1=HDL~AGE+SBP+factor(h.N1AA), 
	formula2=HDL~AGE+SBP, pheno=pheno.dat, haplo=myhaplo, sim=10, 
	predict_variable="h.N1AA")
summary(mymodel)

# example with a subsetting variable, looking at males only
# and modelling a dominant haplotypic effect
mymodel <- haplo.quant(formula1=HDL~AGE+SBP+h.N1AA, formula2=HDL~AGE+SBP, 
	pheno=pheno.dat, haplo=myhaplo, sim=10, effect="dom", 
	sub=expression(SEX==1))



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
