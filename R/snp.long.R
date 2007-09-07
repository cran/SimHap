`snp.long` <- 
function(fixed, random, geno, pheno, cor="corCAR1", form=~1, value=0.2, sub=NULL){
  gennames <- names(geno)
  geno <- data.frame(geno[,-1])
  names(geno) <- gennames[-1]
    
  library(nlme)
  formula1_fixed <- formula(fixed)
  formula_random <- formula(random)
  
  geno_expand <- matrix(0, nrow=nrow(pheno), ncol=ncol(geno))
  geno.mat <- as.matrix(geno)


  de_row <- 1

    geno_expand[1,] <- geno.mat[de_row,]

    last_phen_ID <- pheno[1,1]


    phen_ID <- pheno[,1]
    for(h in 2:nrow(pheno)) {
      current_phen_ID <- phen_ID[h]

      if (current_phen_ID != last_phen_ID) {
        de_row <- de_row + 1
      }

      geno_expand[h,] <- geno.mat[de_row,]
      last_phen_ID <- current_phen_ID
    }
    
    geno_expand <- as.data.frame(geno_expand)
    colnames(geno_expand) <- names(geno)
    dataframe <- as.data.frame(cbind(pheno, geno_expand))
    
    formula1_fixednofactors <- formula1_fixed
    formula1_fixedterms <- attr(terms(formula1_fixednofactors), "term.labels")
    
    if(any(regexpr(":", formula1_fixedterms)!=-1)){
        formula1_fixedterms <- formula1_fixedterms[-which(regexpr(":", formula1_fixedterms)!=-1)]
    }

    if(any(regexpr("factor", formula1_fixedterms)==1)) {
  
        formula1_fixedterms[which(regexpr("factor", formula1_fixedterms)==1)] <- substr(formula1_fixedterms[which(regexpr("factor", formula1_fixedterms)==1)],8,nchar(formula1_fixedterms[which(regexpr("factor", formula1_fixedterms)==1)])-1)

    } 
    #else formula1_fixedterms <- attr(terms(formula1_fixednofactors), "term.labels")
  
  #recreate mydata such that it contains only indivs with complete data
    dataframe <- dataframe[complete.cases(dataframe[formula1_fixedterms]),]

    if(cor=="corCAR1") {
    if(length(sub)==0) {
      fit1.lme <- lme(fixed=fixed, data=dataframe, random=random, correlation=corCAR1(form= form), na.action=na.omit)
      }
    
    else {
      fit1.lme <- eval(substitute(lme(fixed=fixed, data=dataframe, random=random, correlation=corCAR1(value=value, form= form), 
      na.action=na.omit, subset=subset),list(subset=sub)))
      }
    }

    if(cor=="corAR1") {
    if(length(sub)==0) {
      fit1.lme <- lme(fixed=fixed, data=dataframe, random=random, correlation=corAR1(value=value, form= form), na.action=na.omit)
      }
    
    else {
      fit1.lme <- eval(substitute(lme(fixed=fixed, data=dataframe, random=random, correlation=corAR1(value=value, form= form), 
      na.action=na.omit, subset=subset),list(subset=sub)))
      }
    }

 
    if(cor=="corCompSymm") {
    if(length(sub)==0) {
      fit1.lme <- lme(fixed=fixed, data=dataframe, random=random, correlation=corCompSymm(value=value, form=form), na.action=na.omit)
      }
    
    else {
      fit1.lme <- eval(substitute(lme(fixed=fixed, data=dataframe, random=random, correlation=corCompSymm(value=value, form=form), 
      na.action=na.omit, subset=subset),list(subset=sub)))
      }
    }

  anovfull <- as.data.frame(anova(fit1.lme, test="Chisq"))
  anovfull.out <- cbind(anovfull$numDF, signif(anovfull$"p-value", digits=3))

  row.names(anovfull.out) <- row.names(anovfull)
  anovfull.out <- as.data.frame(anovfull.out)
  names(anovfull.out) <- c("DF", "P-Value") 

  lnLbig <- logLik(fit1.lme)
  
  sum.lme <- as.data.frame(summary(fit1.lme)$tTable)
    
  out <- NULL
  out$Coef <- sum.lme$Value
  out$stderror <- sum.lme$Std.Error
  out$pvals <- sum.lme$"p-value"
    
  out <- data.frame(out)
  row.names(out) <- row.names(sum.lme)
  names(out) <- c("Coefficient", "Std.Error", "P.Value")
  
  out.list <- list(results=out, fixed_formula=fixed, random_formula=random, ANOD=anovfull.out, loglik=lnLbig, fit.lme=fit1.lme, AIC=AIC(fit1.lme), corStruct=cor)
  class(out.list) <- "snpLong"
  return(out.list)
  
} 

