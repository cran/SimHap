`epi.surv` <- 
function(formula, pheno, sub=NULL){

  library(survival)
  
  coxph.out <- function(x, digits=max(options()$digits - 4, 3), ...)
    {

      if (!is.null(x$fail)) {
	    cat(" Coxph failed.", x$fail, "\n")
	    return()
	    }
	    
      savedig <- options(digits = digits)
      on.exit(options(savedig))

      coef <- x$coef
      se <- sqrt(diag(x$var))
      if(is.null(coef) | is.null(se))
        stop("Input is not valid")

      if (is.null(x$naive.var)) {
	    tmp <- cbind(coef, exp(coef), se, coef/se,
	       signif(1 - pchisq((coef/ se)^2, 1), digits -1))
	    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	       "se(coef)", "z", "p"))
	    }
      
      else {
	      nse <- sqrt(diag(x$naive.var))
	      tmp <- cbind(coef, exp(coef), nse, se, coef/se,
	        signif(1 - pchisq((coef/se)^2, 1), digits -1))
	      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	        "se(coef)", "robust se", "z", "p"))
	    }
      cat("\n")
      tmp <- as.data.frame(tmp)
      return(tmp)
    }

  if(length(sub)==0) {

      fit1.coxph <- coxph(formula=formula, data=pheno, na.action=na.omit)

      }
    
    else {
      fit1.coxph <- eval(substitute(coxph(formula=formula, data=pheno, na.action=na.omit, subset=subset), list(subset=sub)))
      }

    wald <- formatC(summary(fit1.coxph)$waldtest)
    wald <- as.data.frame(cbind(wald[1], wald[2], wald[3]))
    row.names(wald) <- ""
    names(wald) <- c("Wald test", "df", "P-value")
    out.coxph <- coxph.out(fit1.coxph)

    logLik <- fit1.coxph$loglik
    r.squared <- summary(fit1.coxph)$rsq[1]
    max.rsquared <- summary(fit1.coxph)$rsq[2]
    rsquared.out <- rbind(round(r.squared,digits=3), round(max.rsquared, digits=3))
    rsquared.out <- as.data.frame(rsquared.out)
    names(rsquared.out) <- "R-Squared"
    rownames(rsquared.out) <- c("Model r-squared", "Max r-squared")
    
    out <- NULL
    out$OR <- signif(out.coxph$"exp(coef)", digits=4)
    out$OR.CI <- paste("(",signif(exp(out.coxph$coef - 1.96*(out.coxph$"se(coef)")),digits=3),",",signif(exp(out.coxph$coef + 1.96*(out.coxph$"se(coef)")), digits=3),")", sep="")
    out$pvals <- out.coxph$p
    

    output <- cbind(out$OR, out$OR.CI, out$pvals)
    output <- data.frame(output)
    row.names(output) <- row.names(out.coxph)
    names(output) <- c("HR", "HR.95%CI", "P.Value")
    
    out.list <- list(results=output, formula=formula, fit.coxph=fit1.coxph, Wald=wald, logLik=logLik, rsquared=rsquared.out)
    class(out.list) <- "epiSurv"
    return(out.list)
  

 
}
