\name{SNP.dat}
\alias{SNP.dat}
\docType{data}
\title{Example genotypic data}
\description{
  \code{SNP.dat} is an example genotypic data set containing single nucleotide polymorphism (SNP) genotypes.
}
\usage{data(SNP.dat)}
\format{
  A data frame with 180 observations on the following 5 variables.
  \describe{
    \item{\code{ID}}{patient identifier.}
    \item{\code{SNP_1}}{a SNP vector containing biallelic SNP genotypes.}
    \item{\code{SNP_2}}{a SNP vector containing biallelic SNP genotypes.}
    \item{\code{SNP_3}}{a SNP vector containing biallelic SNP genotypes.}
    \item{\code{SNP_4}}{a SNP vector containing biallelic SNP genotypes.}
  }
}
\details{
  Simulated SNP data for use with \code{pheno.dat}
}

\examples{

data(SNP.dat)

# transforms SNP.dat to an object containing 3 columns
# per SNP - additive, dominant and recessive, where genotypes
# defined in 'baseline' serve as the baseline genotypes
geno.dat <- SNP2Geno(SNP.dat, baseline=c("MM", "11", "GG", "CC"))

data(pheno.dat)
mymodel <- snp.quant(formula1=HDL~AGE+SBP+factor(SNP_1_add),
	formula2=HDL~AGE+SBP, geno=geno.dat, pheno=pheno.dat)
summary(mymodel)

}
\keyword{datasets}
