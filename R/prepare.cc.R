'prepare.cc' <- 
function(geno, pheno, cc.var){

    if(any(is.na(pheno[cc.var]))){
         newgeno <- geno[complete.cases(pheno[cc.var]),]
         newpheno <- pheno[complete.cases(pheno[cc.var]),]
    }
    else {
         newgeno <- geno
         newpheno <- pheno
    }

    out <- list(geno=newgeno, pheno=newpheno)
    return(out)
}
