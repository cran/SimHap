\name{infer.haplos.cc}
\alias{infer.haplos.cc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Infer haplotype configuration independently in cases and controls}
\description{
  \code{infer.haplos.cc} generates a haplotype object to be used in association analysis.
}
\usage{
infer.haplos.cc(geno, pheno, cc.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{a genotype data frame where each SNP is represented by two columns, one for each allele, in the form of \code{haplo.dat}.}
  \item{pheno}{a data frame containing phenotype data with at least two columns - a subject identifier and an indicator of disease status.}
  \item{cc.var}{the column name of the parameter indicating disease status. Must be entered with quotations, e.g. ``DISEASE".}
}
\details{
  \code{cc.var} must be binary, taking only values 0 or 1. 
}
\value{
  \code{infer.haplos.cc} returns a list containing the following items
  \item{hapMat}{a dataframe containing all possible haplotype configurations with their respective likelihoods, for each individual.}
  \item{hap.freq.cases}{haplotype frequencies among cases estimated using the EM algorithm, and the standard errors of these frequencies.}
  \item{hap.freq.controls}{haplotype frequencies among controls estimated using the EM algorithm, and the standard errors of these frequencies.}
  \item{init.freq.cases}{initial haplotype frequencies among cases to be used by other SimHap functions.}
  \item{init.freq.controls}{initial haplotype frequencies among controls to be used by other SimHap functions.}
  
}
\references{
McCaskie, P.A., Carter, K.W, Hazelton, M., Palmer, L.J. (2007) SimHap: A comprehensive modeling framework for epidemiological outcomes and a multiple imputation approach to haplotypic analysis of population-based data, [online] www.genepi.org.au/simhap.

Stram, D.O., Leigh Pearce, C., Bretsky, P., Freedman, M., Hirschhorn, J.N., Altshuler, D., Kolonel, L.N., Henderson, B.E., Thomas, D.C. (2003) Modeling and EM Estimation of Haplotype-Specific Relative Risks from Genotype Data for a Case-Control Study of Unrelated Individuals, \emph{Human Heredity}, 55:179-190.
}
\author{Pamela A. McCaskie}
\note{\code{infer.haplos.cc} is to be used in place of \code{infer.haplos} when haplotypes and haplotype frequencies are to be inferred independently in cases and controls. \code{geno} and \code{pheno} should have individuals in the same order, with the subject identifier column in ascending order.
}
\seealso{\code{\link{infer.haplos}}, \code{\link{prepare.cc}}}
\examples{

data(SNP.dat)

# convert SNP.dat to format required by infer.haplos
haplo.dat <- SNP2Haplo(SNP.dat)

data(pheno.dat)

newdata <- prepare.cc(geno=haplo.dat, pheno=pheno.dat, cc.var="DISEASE")
newhaplo.dat <- newdata$geno
newpheno.dat <- newdata$pheno

# generates haplotype frequencies and haplotype design matrix
myinfer<-infer.haplos.cc(geno=newhaplo.dat, 
	pheno=newpheno.dat, cc.var="DISEASE")

# prints haplotype frequencies among cases
myinfer$hap.freq.cases

# prints haplotype frequencies among controls
myinfer$hap.freq.controls
 
# generated haplo object where haplotypes with a frequency 
# below min.freq are grouped as a category called "rare"
myhaplo<-make.haplo.rare(myinfer,min.freq=0.05) 
mymodel <- haplo.quant(formula1=HDL~AGE+SBP+h.N1AA, formula2=HDL~AGE+SBP, 
	pheno=newpheno.dat, haplo=myhaplo, sim=10)


}
\keyword{manip}

