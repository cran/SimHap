`epi.bin` <- 
function(formula, pheno, sub=NULL){

  fit1.glm <- eval(substitute(glm(formula=formula, data=pheno, family=binomial, subset=subset),list(subset=sub)))
  anov <- anova(fit1.glm, test="Chisq")
  fit.glm <- as.data.frame(summary(fit1.glm)$coefficients)
  loglikelihood <- logLik(fit1.glm)

  anov.out <- cbind(anov$"Resid. Df", anov$Df, formatC(anov$"P(>|Chi|)"))
  anov.out[1,2] <- ""
  anov.out[1,3] <- ""
  row.names(anov.out) <- row.names(anov)
  anov.out <- as.data.frame(anov.out)
  names(anov.out) <- c("Residual DF", "DF", "P-Value")

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
  
  out.list <- list(formula=formula,results=out,fit.glm=fit1.glm, ANOD=anov.out, logLik=loglikelihood, AIC=AIC(fit1.glm))
  class(out.list) <- "epiBin"
  return(out.list)
  
} 

