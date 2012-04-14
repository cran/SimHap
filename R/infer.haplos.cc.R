`infer.haplos.cc` <-
function(geno, pheno, cc.var){

        if(any(is.na(pheno[cc.var]))) {
	stop("\nMissing values in variable ", cc.var, ". Use 'prepare.cc' to clean data before inferring haplotypes.\n", "See 'help(prepare.cc)' for more details.") 
	}
	
	cases <- geno[pheno[cc.var]==1,]
	controls <- geno[pheno[cc.var]==0,]
	infer.cases <- infer.haplos(cases)
	infer.controls <- infer.haplos(controls)

	hapMat <- rbind(infer.cases$hapMat, infer.controls$hapMat)
	hapMat <- hapMat[order(hapMat$ID),]
	hap.freq.cases <- infer.cases$hap.freq
	hap.freq.controls <- infer.controls$hap.freq
	initFreq.cases <- infer.cases$initFreq
	initFreq.controls <- infer.controls$initFreq

	out <- list(hapMat=hapMat, hap.freq.cases=hap.freq.cases, hap.freq.controls=hap.freq.controls, initFreq.cases=initFreq.cases, initFreq.controls=initFreq.controls)
	return(out)
}
