`summary.snpClogit` <- 
function(object, ...)
{
    object.clogit <- object$fit.clogit
    summary.object.clogit <- summary(object.clogit)
######################################################################
    ans <- list(terms=object.clogit$terms, coefficients=object$results, formula=object$formula1, LRT=object$LRT, Wald=object$Wald,rsquared=object$rsquared, residuals=object.clogit$residuals)
    
    class(ans) <- "summary.snpClogit"
    return(ans)
}
