`epi.long` <- 
function(fixed, random, pheno, cor="corCAR1", value=0.2, form=~1, sub=NULL){

  library(nlme)

    if(cor=="corCAR1") {
    if(length(sub)==0) {
      fit1.lme <- lme(fixed=fixed, data=pheno, random=random, correlation=corCAR1(form= form), na.action=na.omit)
      }
    
    else {
      fit1.lme <- eval(substitute(lme(fixed=fixed, data=pheno, random=random, correlation=corCAR1(form= form), na.action=na.omit,      
      subset=subset),list(subset=sub)))
      }
    }

    if(cor=="corAR1") {
    if(length(sub)==0) {
      fit1.lme <- lme(fixed=fixed, data=pheno, random=random, correlation=corAR1(value=value, form= form), na.action=na.omit)
      }
    
    else {
      fit1.lme <- eval(substitute(lme(fixed=fixed, data=pheno, random=random, correlation=corAR1(value=value, form= form), 
      na.action=na.omit, subset=subset),list(subset=sub)))
      }
    }

    if(cor=="corCompSymm") {
    print("corCompSymm")
    if(length(sub)==0) {
      fit1.lme <- lme(fixed=fixed, data=pheno, random=random, correlation=corCompSymm(value=value, form=form), na.action=na.omit)
      }
    
    else {
      fit1.lme <- eval(substitute(lme(fixed=fixed, data=pheno, random=random, correlation=corCompSymm(value=value, form=form), 
      na.action=na.omit, subset=subset),list(subset=sub)))
      }
    }
   
  sum.lme <- as.data.frame(summary(fit1.lme)$tTable)
  loglikelihood <- logLik(fit1.lme)
  out <- NULL
  out$Coef <- sum.lme$Value
  out$stderror <- sum.lme$Std.Error
  out$pvals <- sum.lme$"p-value"
    
  out <- data.frame(out)
  row.names(out) <- row.names(sum.lme)
  names(out) <- c("Coefficient", "Std.Error", "P.Value")
  
  anovfull <- as.data.frame(anova(fit1.lme, test="Chisq"))
  anovfull.out <- cbind(anovfull$numDF, signif(anovfull$"p-value", digits=3))

  row.names(anovfull.out) <- row.names(anovfull)
  anovfull.out <- as.data.frame(anovfull.out)
  names(anovfull.out) <- c("DF", "P-Value") 

  
  out.list <- list(results=out, fixed_formula=fixed, random_formula=random, fit.lme=fit1.lme, ANOD=anovfull.out, logLik=loglikelihood, AIC=AIC(fit1.lme), corStruct=cor)
  class(out.list) <- "epiLong"

  return(out.list)
  
} 
