`snp.cc.match` <- 
function(formula1, formula2, geno, pheno, sub=NULL){
   geno <- geno[,-1]
   library(survival)
   dataframe <- cbind(pheno, geno)

   #define what the output of a coxph should look like
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
       tmp <- as.data.frame(tmp)
       return(tmp)
     }

   
     formula1 <- formula(formula1)
     formula2 <- formula(formula2)
  
     #takes care of variables defined as factors in the formulae so they can be found in the data
     formula1_nofactors <- formula1
  if(any(regexpr("factor", attr(terms(formula1_nofactors), "term.labels"))==1)) {
  
        formula1_terms <- attr(terms(formula1_nofactors), "term.labels")
        formula1_terms[which(regexpr("factor", attr(terms(formula1_nofactors), "term.labels"))==1)] <- substr(attr(terms(formula1_nofactors), "term.labels")[which(regexpr("factor", attr(terms(formula1_nofactors), "term.labels"))==1)],8,nchar(attr(terms(formula1_nofactors), "term.labels")[which(regexpr("factor", attr(terms(formula1_nofactors), "term.labels"))==1)])-1)
	
	

  } 
  
  if(any(regexpr("strata", attr(terms(formula1_nofactors), "term.labels"))==1)){
  
        formula1_terms <- attr(terms(formula1_nofactors), "term.labels")
        formula1_terms[which(regexpr("strata", attr(terms(formula1_nofactors), "term.labels"))==1)] <- substr(attr(terms(formula1_nofactors), "term.labels")[which(regexpr("strata", attr(terms(formula1_nofactors), "term.labels"))==1)],8,nchar(attr(terms(formula1_nofactors), "term.labels")[which(regexpr("strata", attr(terms(formula1_nofactors), "term.labels"))==1)])-1)

  } else formula1_terms <- attr(terms(formula1_nofactors), "term.labels")
  
     
     #dataframe <- as.data.frame(cbind(pheno, geno))

     dataframe <- dataframe[complete.cases(dataframe[formula1_terms]),]

   if(length(sub)==0) {
       fit1.clogit <- clogit(formula=formula1, data=dataframe, na.action=na.omit)
       }

     else {
       fit1.clogit <- eval(substitute(clogit(formula=formula1, data=dataframe, na.action=na.omit, subset=subset), list(subset=sub)))
       }

   if(length(sub)==0) {
       fit2.clogit <- clogit(formula=formula2, data=dataframe, na.action=na.omit)
       }

     else {
       fit2.clogit <- eval(substitute(clogit(formula=formula2, data=dataframe, na.action=na.omit, subset=subset), list(subset=sub)))
       }

   lnLbig <- fit1.clogit$loglik[2]
   lnLsmall <- fit2.clogit$loglik[2]
   fit1.df <- summary(fit1.clogit)$logtest[2]
   fit2.df <- summary(fit2.clogit)$logtest[2]

   lr <- -2*(lnLsmall-lnLbig)
   lr.df <- fit1.df-fit2.df
   lrt <- pchisq(lr,df=lr.df)
   
   anov.out1 <- cbind(round(lnLbig,digits=4), fit1.df, round(lr,digits=4), signif((1-lrt), digits=4))
   anov.out2 <- cbind(round(lnLsmall,digits=4), fit2.df, "", "")
   row.names(anov.out1) <- c("Full model")
   row.names(anov.out2) <- c("Non-genetic")

   LRT <- rbind(anov.out1, anov.out2)
   LRT <- as.data.frame(LRT)
   names(LRT) <- c("logLik", "df", "LR", "P.Value")

   anov <- as.data.frame(anova(fit2.clogit, fit1.clogit, test="Chisq"))
   anovp.row <- anov$"P(>|Chi|)"[2]
   anovdf.row <- anov$Df
   anovloglik <- anov$loglik
   anov.out1 <- cbind(round(anovloglik[1],digits=3), "", "")
   row.names(anov.out1) <- c("1")
   anov.out2 <- cbind(round(anovloglik[2], digits=3), anovdf.row[2], signif(anovp.row, digits=3))
   row.names(anov.out2) <- c("2")
   anov.out <- rbind(anov.out1, anov.out2)
   anov.out <- as.data.frame(anov.out)
   names(anov.out) <- c("loglik", "DF", "P-Value")
 
   fit1.rsquared <- summary(fit1.clogit)$rsq[1]
   fit2.rsquared <- summary(fit2.clogit)$rsq[1]
   max.rsquared <- summary(fit1.clogit)$rsq[2]

   rsquared.out <- rbind(round(fit2.rsquared,digits=3),round(fit1.rsquared,digits=3), round(max.rsquared, digits=3))
   rsquared.out <- as.data.frame(rsquared.out)

   names(rsquared.out) <- "R-Squared"
   row.names(rsquared.out) <- c("Without SNPs", "Including SNPs", "Max R-Squared")

     wald <- formatC(summary(fit1.clogit)$waldtest)
     wald <- as.data.frame(cbind(wald[1], wald[2], wald[3]))
     row.names(wald) <- ""
     names(wald) <- c("Wald test", "df", "P-value")
     out.clogit <- coxph.out(fit1.clogit)

     out <- NULL
     out$OR <- out.clogit$"exp(coef)"
     out$OR.upper <- exp(out.clogit$coef + 1.96*(out.clogit$"se(coef)"))
     out$OR.lower <- exp(out.clogit$coef - 1.96*(out.clogit$"se(coef)"))
     out$pvals <- out.clogit$p

     output <- cbind(out$OR, out$OR.lower, out$OR.upper, out$pvals)
     output <- data.frame(output)
     row.names(output) <- row.names(out.clogit)
     names(output) <- c("OR", "OR.lower", "OR.upper","P.Value")

     out.list <- list(results=output, formula1=formula1, formula2=formula2,LRT=LRT,Wald=wald,logLik=lnLbig[2],fit.clogit=fit1.clogit,fitsub.clogit=fit2.clogit,ANOVA=anov.out,rsquared=rsquared.out)
     class(out.list) <- "snpClogit"
     return(out.list)
  

 }

