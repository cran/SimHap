`summary.hapQuant` <-
function(object, ...)
{
    
######################################################################
    ans <- c(list(formula=object$formula1, coefficients=object$results, empiricalResults=object$empiricalResults, summary.coefs=object$summary.coefs, rsquared=object$rsquared,LRT=object$LRT, AIC=object$aic, predicted=object$predicted, effect=object$effect))
    
    class(ans) <- "summary.hapQuant"
    return(ans)
}

