\name{summary.hapSurv}
\alias{summary.hapSurv}
\alias{print.summary.hapSurv}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize haplotype analysis models for survival data}
\description{
  Summary method for objects of class \code{hapSurv}
}
\usage{
\method{summary}{hapSurv}(object, ...)
\method{print}{summary.hapSurv}(x, digits = max(3, getOption("digits") - 3), 
	signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{hapSurv}, a result of a call to \code{haplo.surv}.}
  \item{x}{an object of class \code{summary.hapSurv}, the result of a call to \preformatted{summary.hapSurv.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{

  \code{summary.hapSurv} returns an object of class \code{summary.hapSurv}, a list with components
  
  \item{formula}{formula1 passed into \code{haplo.surv}.}
  \item{coefficients}{summarized results from fitted model, including hazard ratios and p-values.}
  \item{empiricalResults}{a list containing the hazard ratios, confidence intervals and p-values calculated at each simulation of \code{haplo.surv}.}
  \item{rsquared}{r-squared values for the model including haplotypes and the model excluding haplotypes.}
  \item{LRT}{a likelihood ratio test, testing for significant improvement of the model when haplotypic parameters are included.}
  \item{Wald}{The Wald test for overall significance of the fitted model.}
  \item{effect}{the haplotypic effect modelled: `ADDITIVE', `DOMINANT' or `RECESSIVE'.}
  

}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A. McCaskie}

\seealso{\code{\link{haplo.surv}}}
\examples{

data(SNPsurv.dat)

# convert SNP.dat to format required by infer.haplos
survHaplo.dat <- SNP2Haplo(SNPsurv.dat)

data(survPheno.dat)

# generate haplotype frequencies and haplotype design matrix
myinfer<-infer.haplos(survHaplo.dat) 

# print haplotype frequencies generated by infer.haplos
myinfer$hap.freq 
myhaplo<-make.haplo.rare(myinfer,min.freq=0.05)
mymodel <- haplo.surv(formula1=Surv(time, status)~age+h.V1AA, 
	formula2=Surv(time, status)~age, haplo=myhaplo, 
	pheno=survPheno.dat, sim=10)
summary(mymodel)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

