\name{SNPsurv.dat}
\alias{SNPsurv.dat}
\docType{data}
\title{Example genotypic data for survival analysis}
\description{
  \code{SNPsurv.dat} is an example genotypic data set containing single nucleotide polymorphism (SNP) genotypes for use with \code{snp.surv} and \code{haplo.surv}
}
\usage{data(SNPsurv.dat)}
\format{
  A data frame with 38 observations on the following 4 variables.
  \describe{
    \item{\code{ID}}{patient identifier.}
    \item{\code{SNP_1}}{a SNP vector containing biallelic SNP genotypes.}
    \item{\code{SNP_2}}{a SNP vector containing biallelic SNP genotypes.}
    \item{\code{SNP_3}}{a SNP vector containing biallelic SNP genotypes.}
  }
}
\details{
  Simulated SNP data for use with \code{survPheno.dat}
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
