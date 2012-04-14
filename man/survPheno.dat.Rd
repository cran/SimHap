\name{survPheno.dat}
\alias{survPheno.dat}
\docType{data}
\title{Example phenotypic survival data}
\description{
  Data on the recurrence times to infection, at the point of insertion of the catheter, for kidney patients using portable dialysis equipment. Catheters may be removed for reasons other than infection, in which case the observation is censored.  Each patient has exactly 2 observations.
  \code{pheno.dat} is an example phenotypic data set containing biological measures to be used by \code{snp.surv}, \code{haplo.surv}.
}
\usage{data(survPheno.dat)}
\format{
  A data frame with 76 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{patient identifier.}
    \item{\code{time}}{time.}
    \item{\code{status}}{event status.}
    \item{\code{age}}{in years.}
    \item{\code{sex}}{1=male, 2=female.}
    \item{\code{disease}}{disease type (0=GN, 1=AN, 2=PKD, 3=Other).}
  }
}
\details{
  This data is the same as the \code{kidney} data from the \code{survival} package, with the \code{frail} parameter eliminated.
}
\source{
  McGilchrist and Aisbett, Biometrics 47, 461-66, 1991
}

\examples{

data(SNPsurv.dat)

# transforms SNPlong.dat to an object containing 3 columns
# per SNP - additive, dominant and recessive, where genotypes
# defined in 'baseline' serve as the baseline genotypes
survGeno.dat <- SNP2Geno(SNPsurv.dat, baseline=c("V2V2", "GG", "CC"))

data(survPheno.dat)

mymodel <- snp.surv(formula1=Surv(time, status)~age+SNP_1_add, 
	formula2=Surv(time, status)~age, pheno=survPheno.dat, 
	geno=survGeno.dat)
summary(mymodel)

}
\keyword{datasets}
