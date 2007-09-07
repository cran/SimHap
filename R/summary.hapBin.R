`summary.hapBin` <-
function(object, ...)
{
    
######################################################################
    ans <- c(list(formula=object$formula1, coefficients=object$results,LRT=object$LRT, AIC=object$aic, empiricalResults=object$empiricalResults, summary.coefs=object$summary.coefs,effect=object$effect))
    
    class(ans) <- "summary.hapBin"
    return(ans)
}

