`snp.surv` <- 
function(formula1, formula2, geno, pheno, sub=NULL){
   gennames <- names(geno)
   geno <- data.frame(geno[,-1])
   names(geno) <- gennames[-1]
    
   library(survival)

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
     
     formula1 <- formula(formula1)
     formula2 <- formula(formula2)
  
     #takes care of variables defined as factors in the formulae so they can be found in the data
     formula1_nofactors <- formula1
     formula1_terms <- attr(terms(formula1_nofactors), "term.labels")

     if(any(regexpr(":", formula1_terms)!=-1)){
        formula1_terms <- formula1_terms[-which(regexpr(":", formula1_terms)!=-1)]
     }
     
     if(any(regexpr("factor", formula1_terms)==1)) {
        formula1_terms[which(regexpr("factor", formula1_terms)==1)] <- substr(formula1_terms[which(regexpr("factor", formula1_terms)==1)],8,nchar(formula1_terms[which(regexpr("factor", formula1_terms)==1)])-1)
     } 
   #else formula1_terms <- attr(terms(formula1_nofactors), "term.labels")

     geno_expand <- as.data.frame(geno_expand)
     colnames(geno_expand) <- names(geno)
   
     dataframe <- as.data.frame(cbind(pheno, geno_expand))

     dataframe <- dataframe[complete.cases(dataframe[formula1_terms]),]

   if(length(sub)==0) {
       fit1.coxph <- coxph(formula=formula1, data=dataframe, na.action=na.omit)

       }

     else {
       fit1.coxph <- eval(substitute(coxph(formula=formula1, data=dataframe, na.action=na.omit, subset=subset), list(subset=sub)))
       }

   if(length(sub)==0) {
       fit2.coxph <- coxph(formula=formula2, data=dataframe, na.action=na.omit)
       }

     else {
       fit2.coxph <- eval(substitute(coxph(formula=formula2, data=dataframe, na.action=na.omit, subset=subset), list(subset=sub)))
       }

   lnLbig <- fit1.coxph$loglik[2]
   lnLsmall <- fit2.coxph$loglik[2]
   fit1.df <- summary(fit1.coxph)$logtest[2]
   fit2.df <- summary(fit2.coxph)$logtest[2]
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

   anov <- as.data.frame(anova(fit2.coxph, fit1.coxph, test="Chisq"))
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
   fit1.rsquared <- summary(fit1.coxph)$rsq[1]
   fit2.rsquared <- summary(fit2.coxph)$rsq[1]
   max.rsquared <- summary(fit1.coxph)$rsq[2]

   rsquared.out <- rbind(round(fit2.rsquared,digits=3),round(fit1.rsquared,digits=3), round(max.rsquared, digits=3))
   rsquared.out <- as.data.frame(rsquared.out)

   names(rsquared.out) <- "R-Squared"
   row.names(rsquared.out) <- c("Without SNPs", "Including SNPs", "Max R-Squared")

     wald <- formatC(summary(fit1.coxph)$waldtest)
     wald <- as.data.frame(cbind(wald[1], wald[2], wald[3]))
     row.names(wald) <- ""
     names(wald) <- c("Wald test", "df", "P-value")
     out.coxph <- coxph.out(fit1.coxph)

     out <- NULL
     out$OR <- out.coxph$"exp(coef)"
     out$OR.upper <- exp(out.coxph$coef + 1.96*(out.coxph$"se(coef)"))
     out$OR.lower <- exp(out.coxph$coef - 1.96*(out.coxph$"se(coef)"))
     
     out$pvals <- out.coxph$p

     output <- cbind(out$OR, out$OR.lower, out$OR.upper, out$pvals)
     output <- data.frame(output)
     row.names(output) <- row.names(out.coxph)
     names(output) <- c("HR", "HR.lower", "HR.upper", "P.Value")

     out.list <- list(results=output, formula1=formula1, formula2=formula2,LRT=LRT,Wald=wald,logLik=lnLbig[2],fit.coxph=fit1.coxph,fitsub.coxph=fit2.coxph,ANOVA=anov.out,rsquared=rsquared.out)
     class(out.list) <- "snpSurv"
     return(out.list)
  

 }