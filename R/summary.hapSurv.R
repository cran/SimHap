`summary.hapSurv` <- 
function(object, ...)
{
    
######################################################################
    ans <- c(list(formula=object$formula1, coefficients=object$results, empiricalResults=object$empiricalResults, summary.coefs=object$summary.coefs, rsquared=object$rsquared, LRT=object$LRT, Wald=object$Wald, effect=object$effect))
    
    class(ans) <- "summary.hapSurv"
    return(ans)
}
