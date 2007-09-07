`summary.hapClogit` <-
function(object, ...)
{
    
######################################################################
    ans <- c(list(formula=object$formula1, coefficients=object$results, empiricalResults=object$empiricalResults, summary.coefs=object$summary.coefs, LRT=object$LRT, Wald=object$Wald, rsquared=object$rsquared, effect=object$effect))
    
    class(ans) <- "summary.hapClogit"
    return(ans)
}

