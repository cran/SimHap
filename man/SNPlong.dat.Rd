\name{SNPlong.dat}
\alias{SNPlong.dat}
\docType{data}
\title{Example genotypic data for longitudinal analysis}
\description{
  \code{SNPlong.dat} is an example genotypic data set containing single nucleotide polymorphism (SNP) genotypes for use with \code{snp.long} and \code{haplo.long}
}
\usage{data(SNPlong.dat)}
\format{
  A data frame with 99 observations on the following 4 variables.
  \describe{
    \item{\code{ID}}{patient identifier.}
    \item{\code{SNP_1}}{a SNP vector containing biallelic SNP genotypes.}
    \item{\code{SNP_2}}{a SNP vector containing biallelic SNP genotypes.}
    \item{\code{SNP_3}}{a SNP vector containing biallelic SNP genotypes.}
  }
}
\details{
  Simulated SNP data for use with \code{longPheno.dat}
}

\examples{

data(SNPlong.dat)

# transforms SNPlong.dat to an object containing 3 columns
# per SNP - additive, dominant and recessive, where genotypes
# defined in 'baseline' serve as the baseline genotypes
longGeno.dat <- SNP2Geno(SNPlong.dat, baseline=c("AA", "GG", "V2V2"))

data(longPheno.dat)
mymodel <- snp.long(fixed=fev1f~SNP_1_add, random=~1|id, 
	geno=longGeno.dat, pheno=longPheno.dat, form=~year|id)
summary(mymodel)

}
\keyword{datasets}
