`epi.cc.match` <- 
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

      fit1.clogit <- clogit(formula=formula, data=pheno, na.action=na.omit)

      }
    
    else {
      fit1.clogit <- eval(substitute(clogit(formula=formula, data=pheno, na.action=na.omit, subset=subset), list(subset=sub)))
      }

    wald <- formatC(summary(fit1.clogit)$waldtest)
    wald <- as.data.frame(cbind(wald[1], wald[2], wald[3]))
    row.names(wald) <- ""
    names(wald) <- c("Wald test", "df", "P-value")
    out.clogit <- coxph.out(fit1.clogit)

    logLik <- fit1.clogit$loglik
    r.squared <- summary(fit1.clogit)$rsq[1]
    max.rsquared <- summary(fit1.clogit)$rsq[2]
    rsquared.out <- rbind(round(r.squared,digits=3), round(max.rsquared, digits=3))
    rsquared.out <- as.data.frame(rsquared.out)
    names(rsquared.out) <- "R-Squared"
    rownames(rsquared.out) <- c("Model r-squared", "Max r-squared")
    
    out <- NULL
    out$OR <- out.clogit$"exp(coef)"
    out$OR.upper <- exp(out.clogit$coef + 1.96*(out.clogit$"se(coef)"))
    out$OR.lower <- exp(out.clogit$coef - 1.96*(out.clogit$"se(coef)"))
    out$pvals <- out.clogit$p
    

    output <- cbind(out$OR, out$OR.lower, out$OR.upper, out$pvals)
    output <- data.frame(output)
    row.names(output) <- row.names(out.clogit)
    names(output) <- c("OR", "OR.lower", "OR.upper","P.Value")
    
    out.list <- list(results=output, formula=formula, fit.clogit=fit1.clogit, Wald=wald, logLik=logLik, rsquared=rsquared.out)
    class(out.list) <- "epiClogit"
    return(out.list)
  

 
}
