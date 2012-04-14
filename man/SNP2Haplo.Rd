\name{SNP2Haplo}
\alias{SNP2Haplo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert biallelic SNP data to haplo format}
\description{
  Converts biallelic SNP data to the format required by the suite of haplotype inference and analysis functions
}
\usage{
SNP2Haplo(geno)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{A data frame containing biallelic SNP data. The first column should contain a subject identifier. The remaining columns should each represent a biallelic SNP, containing two alleles per entry with no separator.}
}

\value{
  \code{SNP2Haplo} returns a data frame whose first column is a subject identifier, and where the remaining columns represent the genotype data. Two columns are returned per SNP, each containing a single allele. 
}

\author{Pamela A. McCaskie}

\seealso{\code{\link{SNP2Geno}}}
\examples{

data(SNP.dat)
haplo.dat <- SNP2Haplo(SNP.dat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
