`summary.epiClogit` <- 
function(object, ...)
{
    object.clogit <- object$fit.clogit
    summary.object.clogit <- summary(object.clogit)

    ans <- list(terms=object.clogit$terms, coefficients=object$results, formula=object$formula, Wald=object$Wald, logLik=object$logLik,rsquared=object$rsquared, residuals=object.clogit$residuals)
    
    class(ans) <- "summary.epiClogit"
    return(ans)
}
