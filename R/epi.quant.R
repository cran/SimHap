`epi.quant` <- 
function(formula, pheno, sub=NULL){


  fit1.lm <- eval(substitute(lm(formula=formula, data=pheno, subset=subset),list(subset=sub)))
  
  anov <- as.data.frame(anova(fit1.lm, test="Chisq"))
  fit.lm <- as.data.frame(summary(fit1.lm)$coefficients)
  loglikelihood <- logLik(fit1.lm)
 
  anovfull <- as.data.frame(anova(fit1.lm, test="Chisq"))
  anovfull.out <- cbind(anovfull$Df, formatC(anovfull$"Pr(>F)"))
  #anovfull.out[1,2] <- ""
  #anovfull.out[1,3] <- ""

  row.names(anovfull.out) <- row.names(anovfull)
  anovfull.out <- as.data.frame(anovfull.out)
  names(anovfull.out) <- c("DF", "P-Value")
  
  out <- NULL
  out$coef <- formatC(fit.lm$Estimate)
  out$se <- formatC(fit.lm$"Std. Error")
  out$pval <- formatC(as.vector(t(fit.lm[,ncol(fit.lm)])))
  
  out <- data.frame(out)
  row.names(out) <- row.names(fit.lm)
  names(out) <- c("Coefficient", "Std.Error", "P.Value")
  
  "%w/o%" <- function(x,y) x[!x %in% y]
  invars <- names(fit1.lm$coef)
  check <- invars %w/o% row.names(out)
  if(length(check) != 0) cat(c(check, "removed due to singularities"), "\n")
  
  out.list <- list(formula=formula,results=out,fit.lm=fit1.lm, ANOD=anovfull.out, logLik=loglikelihood, AIC=AIC(fit1.lm)) 
  class(out.list) <- "epiQuant"
  return(out.list)
  }
 