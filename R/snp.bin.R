`snp.bin` <-
function(formula1, formula2, geno, pheno, sub=NULL){
  gennames <- names(geno)
  geno <- data.frame(geno[,-1])
  names(geno) <- gennames[-1]
  mydata <- cbind(pheno, geno)
  
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
 
  #recreate mydata such that it contains only indivs with complete data
  mydata <- mydata[complete.cases(mydata[formula1_terms]),] 

  fit1.glm <- eval(substitute(glm(formula=formula1, data=mydata, family=binomial, subset=subset),list(subset=sub)))
  fit2.glm <- eval(substitute(glm(formula=formula2, data=mydata, family=binomial, subset=subset),list(subset=sub)))
  
  anovfull <- as.data.frame(anova(fit1.glm, test="Chisq"))
  anovfull.out <- cbind(anovfull$"Resid. Df",anovfull$Df, formatC(anovfull$"P(>|Chi|)"))
  anovfull.out[1,2] <- ""
  anovfull.out[1,3] <- ""

  row.names(anovfull.out) <- row.names(anovfull)
  anovfull.out <- as.data.frame(anovfull.out)
  names(anovfull.out) <- c("Residual DF", "DF", "P-Value")

  lnLbig <- logLik(fit1.glm)
  lnLsmall <- logLik(fit2.glm)

  lr <- -2*(lnLsmall[1]-lnLbig[1])
  lr.df <- attr(lnLbig, "df")-attr(lnLsmall,"df")
  lrt <- pchisq(lr,df=lr.df)

  fit.glm <- as.data.frame(summary(fit1.glm)$coefficients)
  
  anov.out1 <- cbind(round(lnLbig[1],digits=4), attr(lnLbig, "df"), round(lr,digits=4), signif((1-lrt), digits=4))
  anov.out2 <- cbind(round(lnLsmall[1],digits=4), attr(lnLsmall, "df"), "", "")
  row.names(anov.out1) <- c("Full model")
  row.names(anov.out2) <- c("Non-genetic")

  LRT <- rbind(anov.out1, anov.out2)
  LRT <- as.data.frame(LRT)
  names(LRT) <- c("logLik", "df", "LR", "P-Value")
  
  out <- NULL
  out$OR <- exp(fit.glm$Estimate)
  out$ORlower <- exp(fit.glm$Estimate - 1.96*(fit.glm$"Std. Error"))
  out$ORupper <- exp(fit.glm$Estimate + 1.96*(fit.glm$"Std. Error"))
  out$pval <- as.vector(t(fit.glm[,ncol(fit.glm)]))
  
  out <- data.frame(out)
  row.names(out) <- row.names(fit.glm)
  names(out) <- c("OR", "ORlower", "ORupper","P.Value")
  
  
  "%w/o%" <- function(x,y) x[!x %in% y]
  invars <- names(fit1.glm$coef)
  check <- invars %w/o% row.names(out)
  if(length(check) != 0) cat(c(check, "removed due to singularities"), "\n")

  out.list <- list(results=out, formula1=formula1, formula2=formula2,LRT=LRT,ANOD=anovfull.out,logLik=lnLbig,fit.glm=fit1.glm,fitsub.glm=fit2.glm,AIC=AIC(fit1.glm))
  class(out.list) <- "snpBin"
  return(out.list)
  
}

