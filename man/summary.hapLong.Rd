\name{summary.hapLong}
\alias{summary.hapLong}
\alias{print.summary.hapLong}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing single SNP analysis models for longitudinal data}
\description{
Summary method for objects of class \code{snpLong}
}
\usage{
\method{summary}{hapLong}(object, ...)
\method{print}{summary.hapLong}(x, digits = max(3, getOption("digits") - 3),
	signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{hapLong}, a result of a call to \code{haplo.long}.}
  \item{x}{an object of class \code{summary.hapLong}, the result of a call to \preformatted{summary.hapLong.}}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, ``significance stars" are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.hapLong} returns an object of class summary.hapLong. A list with components
  \item{fixed_formula}{fixed effects formula used in \code{haplo.long}.}
  \item{random_formula}{random effects formula used in \code{haplo.long}.}
  \item{coefficients}{summarized results from fitted model, including coefficients, standard errors and p-values.}
  \item{empiricalResults}{a list containing the coefficients, standard errors and p-values calculated at each simulation of \code{haplo.long}.}
  \item{AIC}{Akaike Information Criterion for the model fitted in \code{haplo.long}.}
  \item{corStruct}{the correlation structure specified in \code{haplo.long}.}
  \item{effect}{the haplotypic effect modelled: `ADDITIVE', `DOMINANT' or `RECESSIVE'.}
}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A. McCaskie}

\seealso{\code{\link{haplo.long}}}
\examples{

data(SNPlong.dat)

# convert SNP.dat to format required by infer.haplos
longHaplo.dat <- SNP2Haplo(SNPlong.dat)

data(longPheno.dat)

myinfer<-infer.haplos(longHaplo.dat)
myinfer$hap.freq

myhaplo<-make.haplo.rare(myinfer,min.freq=0.05)
mymodel <- haplo.long(fixed=fev1f~h.ACV2, random=~1|ID,
	pheno=longPheno.dat, haplo=myhaplo, cor="corCAR1", form=~year|ID, sim=10)
summary(mymodel)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
