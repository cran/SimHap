\name{longPheno.dat}
\alias{longPheno.dat}
\docType{data}
\title{Example phenotypic, longitudinal data}
\description{
  \code{longPheno.dat} is an example phenotypic data set containing biological measures to be used by \code{snp.long} and \code{haplo.long}
}
\usage{data(longPheno.dat)}
\format{
  A data frame with 601 observations on the following 12 variables.
  \describe{
    \item{\code{id}}{patient identifier.}
    \item{\code{year}}{year of survey.}
    \item{\code{time}}{time point in years, where 1966=1.}
    \item{\code{sex}}{1=male, 0=female.}
    \item{\code{age}}{age in years.}
    \item{\code{height}}{height in metres.}
    \item{\code{weight}}{weight in kilograms.}
    \item{\code{bmi}}{body-mass index.}
    \item{\code{fev1f}}{forced expired volume in the first second - measure of lung function.}
  }
}
\details{
  \code{longPheno.dat} was simulated to take the format of the Busselton Health Survey from Western Australia. 
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
