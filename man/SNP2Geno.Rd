\name{SNP2Geno}
\alias{SNP2Geno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert biallelic SNP data to genotype count data}
\description{
  Converts biallelic SNP data to the format required by SNP analysis functions
}
\usage{
SNP2Geno(geno, baseline)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{A data frame containing biallelic SNP data. The first column should contain a subject identifier. The remaining columns should each represent a biallelic SNP, containing two alleles per entry with no separator.}
  \item{baseline}{A vector of character strings representing the 'wildtype' genotype.}
}
\details{
  \code{SNP2Geno} uses \code{baseline} to recode biallelic SNP data into a series of columns representing genotype counts under additive, dominant and recessive models. The genotype specified as the wildtype will always be recoded as 0. Under an additive model, the heterozygote is recoded as 1 and the minor homozygote is recoded as 2. Under a dominant model, both the heterozygote and minor homozygote are recoded as 1, and under a recessive model, the heterozygote is recoded as 0 and the minor homozygote is coded as 1.
}
\value{
  \code{SNP2Geno} returns a data frame whose first column is a subject identifier. The remaining columns comprise 3 per SNP where the first is recoded under an additive model, the second is recoded under a dominant model and the third is recoded under a recessive model.
}

\author{Pamela A. McCaskie}

\seealso{\code{\link{SNP2Haplo}}}
\examples{

data(SNP.dat)
geno.dat <- SNP2Geno(SNP.dat, baseline=c("MM", "11", "GG", "CC"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

