`summary.epiSurv` <- 
function(object, ...)
{
    object.coxph <- object$fit.coxph
    summary.object.coxph <- summary(object.coxph)

    ans <- list(terms=object.coxph$terms, coefficients=object$results, formula=object$formula, Wald=object$Wald, logLik=object$logLik,rsquared=object$rsquared, residuals=object.coxph$residuals)
    
    class(ans) <- "summary.epiSurv"
    return(ans)
}
