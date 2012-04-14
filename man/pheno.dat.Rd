\name{pheno.dat}
\alias{pheno.dat}
\docType{data}
\title{Example phenotypic data}
\description{
  \code{pheno.dat} is an example phenotypic data set containing biological measures to be used by \code{snp.quant}, \code{haplo.quant}, \code{snp.bin} and \code{haplo.bin}
}
\usage{data(pheno.dat)}
\format{
  A data frame with 180 observations on the following 16 variables.
  \describe{
    \item{\code{ID}}{patient identifiers.}
    \item{\code{SEX}}{1=male, 0=female.}
    \item{\code{AGE}}{age in years.}
    \item{\code{SBP}}{systolic blood pressure (mmHg).}
    \item{\code{DBP}}{diastolic blood pressure (mmHg).}
    \item{\code{BMI}}{body-mass index.}
    \item{\code{WHR}}{waist-hip ratio.}
    \item{\code{HDL}}{plasma high density lipoprotein (mmol/L).}
    \item{\code{LDL}}{plasma low density lipoprotein (mmol/L).}
    \item{\code{DIABETES}}{a binary indicator of history of type 2 diabetes.}
    \item{\code{FH_IHD}}{a binary indicator of family history of ischaemic heart disease.}
    \item{\code{PLAQUE}}{a binary indicator of the presence of 1 or more carotid plaques.}
    \item{\code{SMOKE}}{a binary indicator of smoking history (0=never smoked, 1=ever smoked).}
    \item{\code{PY}}{pack-years of smoking.}
    \item{\code{DISEASE}}{a binary indicator of ischaemic heart disease.}
    \item{\code{STRAT}}{a matching variable indicating the pairs of matched cases and controls.}
  }
}
\details{
    \code{pheno.dat} is a simulated data set of coronary heart disease related phenotypes
}

\examples{

data(SNP.dat)

# convert SNP.dat to format required by snp.quant
geno.dat <- SNP2Geno(SNP.dat, baseline=c("MM", "11", "GG", "CC"))

data(pheno.dat)

mymodel <- snp.quant(formula1=LDL~AGE+SBP+factor(SNP_1_add),
	formula2=HDL~AGE+SBP, geno=geno.dat, pheno=pheno.dat)
summary(mymodel)

# example with a subsetting variable, looking at males only
mymodel <- snp.quant(formula1=LDL~AGE+SBP+factor(SNP_1_add),
	formula2=HDL~AGE+SBP, geno=geno.dat, pheno=pheno.dat, 
	sub=expression(SEX==1))
 
}
\keyword{datasets}
