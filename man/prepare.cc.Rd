\name{prepare.cc}
\alias{prepare.cc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prepare case-control data for inferring haplotypes}
\description{
\code{prepare.cc} prepares case-control data when there may be missing values in the `case status' variable. This eliminates problems when using \code{infer.haplos.cc}.
}
\usage{
prepare.cc(geno, pheno, cc.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{a genotype data frame where each SNP is represented by two columns, one for each allele, in the form of \code{haplo.dat}.}
  \item{pheno}{a data frame containing phenotype data with at least two columns - a subject identifier and an indicator of disease status.}
  \item{cc.var}{the column name of the parameter indicating disease status. Must be entered with quotations, e.g. ``DISEASE".}
}
\details{
  \code{prepare.cc} searches for missing values in \code{cc.var} and reduces \code{geno} and \code{pheno} to include only those individuals with known disease status. These \code{geno} and \code{pheno} objects can then be passed into \code{infer.haplos.cc}.
}
\value{
  
  \item{geno}{a genotype data frame where each SNP is represented by two columns, one for each allele, in the form of \code{haplo.dat}. Individuals with unknown disease status are removed.}
  \item{pheno}{a data frame containing phenotype data with at least two columns - a subject identifier and an indicator of disease status. Individuals with unknown disease status are removed.}

}
\references{
McCaskie, P.A., Carter, K.W. Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple-imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.
}
\author{Pamela A. McCaskie}

\seealso{\code{\link{infer.haplos.cc}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by infer.haplos.cc
haplo.dat <- SNP2Haplo(SNP.dat)
data(pheno.dat)

# not run: will return an error due to missing data in variable 'DISEASE'
# myinfer<-infer.haplos.cc(geno=haplo.dat, pheno=pheno.dat, 
# 	cc.var="DISEASE") 

newdata <- prepare.cc(geno=haplo.dat, pheno=pheno.dat, cc.var="DISEASE")
newhaplo.dat <- newdata$geno
newpheno.dat <- newdata$pheno
myinfer<-infer.haplos.cc(geno=newhaplo.dat, pheno=newpheno.dat, 
	cc.var="DISEASE")

# prints haplotype frequencies among cases
myinfer$hap.freq.cases

# prints haplotype frequencies among controls
myinfer$hap.freq.controls 

# generated haplo object where haplotypes with a frequency 
# below min.freq are grouped as a category called "rare"
myhaplo<-make.haplo.rare(myinfer,min.freq=0.05) 
mymodel <- haplo.bin(formula1=DISEASE~AGE+SBP+h.N1AA, 
	formula2=DISEASE~AGE+SBP, pheno=newpheno.dat, haplo=myhaplo, 
	sim=10)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

