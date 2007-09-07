`summary.snpSurv` <- 
function(object, ...)
{
    object.coxph <- object$fit.coxph
    summary.object.coxph <- summary(object.coxph)
######################################################################
    ans <- list(terms=object.coxph$terms, coefficients=object$results, formula=object$formula1, LRT=object$LRT, Wald=object$Wald,rsquared=object$rsquared, residuals=object.coxph$residuals)
    
    class(ans) <- "summary.snpSurv"
    return(ans)
}
