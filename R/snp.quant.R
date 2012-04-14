`snp.quant` <-
function(formula1, formula2, geno, pheno, sub=NULL, predict_variable=NULL){

  if(!identical(as.character(unique(pheno$ID)), as.character(unique(geno$ID)))) stop("Phenotype data and Genotype data are not in the same order.")

  gennames <- names(geno)
  geno <- data.frame(geno[,-1])
  names(geno) <- gennames[-1]
  mydata <- cbind(pheno, geno)


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

  mydata <- mydata[complete.cases(mydata[formula1_terms]),]

  fit1.lm <- eval(substitute(lm(formula=formula1, data=mydata, subset=subset),list(subset=sub)))
  fit2.lm <- eval(substitute(lm(formula=formula2, data=mydata, subset=subset),list(subset=sub)))

  if(is.null(predict_variable)){
          predicted.values <- "No variable was chosen to generate marginal means"
          }

  else if(!is.null(predict_variable)) {

          nls <- length(levs <- levels(as.factor(mydata[, predict_variable])))
	  tmp1 <- colMeans(mydata, na.rm=TRUE)
          tmp2 <- data.frame(matrix(rep(tmp1, nls), nrow=nls, byrow=TRUE))
	  names(tmp2) <- names(tmp1)
          tmp2[, predict_variable] <- as.factor(levs)
	  mydata[, predict_variable] <- as.factor(mydata[, predict_variable])
	  formula.pred <- formula(paste(rownames(attr(terms(formula1),"factors"))[1],"~", paste(formula1_terms[1:length(formula1_terms)],collapse="+")))
  	  fit.pred <- eval(substitute(lm(formula=formula.pred, data=mydata, subset=subset),list(subset=sub)))

          if(attr(fit1.lm$terms,"dataClasses")[grep(predict_variable, names(attr(fit1.lm$terms,"dataClasses")))]=="factor") {
	  	    #attr(fit1.lm$terms,"dataClasses")[-grep(predict_variable, names(attr(fit1.lm$terms,"dataClasses")))] <- "numeric"
		    predicted.values <- predict(fit.pred, newdata=tmp2, se.fit=T)
                    predicted.values <- cbind(predicted.values$fit,predicted.values$se.fit)
		    if(nrow(predicted.values)==3) row.names(predicted.values) <- c("Wildtype","Heterozygote","Mutation Homozygote")
		    if(nrow(predicted.values)==2) row.names(predicted.values) <- c("Wildtype","Mutation")
                    colnames(predicted.values) <- c("Predicted Value", "Std. Error")
                    } else predicted.values <- "SNP variable was fitted with class 'numeric' but class 'factor' was required"
  } else predicted.values <- NULL

  anovfull <- as.data.frame(anova(fit1.lm, test="Chisq"))
  anovfull.out <- cbind(anovfull$Df, formatC(anovfull$"Pr(>F)"))

  row.names(anovfull.out) <- row.names(anovfull)
  anovfull.out <- as.data.frame(anovfull.out)
  names(anovfull.out) <- c("DF", "P-Value")

  lnLbig <- logLik(fit1.lm)
  lnLsmall <- logLik(fit2.lm)

  lr <- -2*(lnLsmall[1]-lnLbig[1])
  lr.df <- attr(lnLbig, "df")-attr(lnLsmall,"df")
  lrt <- pchisq(lr,df=lr.df)

  fit.lm <- as.data.frame(summary(fit1.lm)$coefficients)

  anov.out1 <- cbind(round(lnLbig[1],digits=4), attr(lnLbig, "df"), round(lr,digits=4), signif((1-lrt), digits=4))
  anov.out2 <- cbind(round(lnLsmall[1],digits=4), attr(lnLsmall, "df"), "", "")
  row.names(anov.out1) <- c("Full model")
  row.names(anov.out2) <- c("Non-genetic")

  LRT <- rbind(anov.out1, anov.out2)
  LRT <- as.data.frame(LRT)
  names(LRT) <- c("logLik", "df", "LR", "P.Value")
##############
  rsquared.out <- rbind(round((summary(fit2.lm)$adj.r.squared),digits=3),round((summary(fit1.lm)$adj.r.squared),digits=3))
  rsquared.out <- as.data.frame(rsquared.out)
  names(rsquared.out) <- "Adjusted R-Squared"
  row.names(rsquared.out) <- c("Without SNP(s)", "Including SNP(s)")
##############
  out <- NULL
  out$coef <- fit.lm$Estimate
  out$se <- fit.lm$"Std. Error"
  out$pval <- as.vector(t(fit.lm[,ncol(fit.lm)]))

  out <- data.frame(out)
  row.names(out) <- row.names(fit.lm)
  names(out) <- c("Coefficient", "Std.Error", "P.Value")

  "%w/o%" <- function(x,y) x[!x %in% y]
  invars <- names(fit1.lm$coef)
  check <- invars %w/o% row.names(out)
  if(length(check) != 0) cat(c(check, "removed due to singularities"), "\n")

  if(is.numeric(predicted.values)) predicted.values <- round(predicted.values, digits=max(3, getOption("digits") - 3))

  out.list <- list(results=out, formula1=formula1, formula2=formula2,LRT=LRT,ANOD=anovfull.out,loglik=lnLbig,fit.lm=fit1.lm,fitsub.lm=fit2.lm,rsquared=rsquared.out,predicted.values=predicted.values, AIC=AIC(fit1.lm))
  class(out.list) <- "snpQuant"
  return(out.list)

}
