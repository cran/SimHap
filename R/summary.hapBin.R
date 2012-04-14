`summary.hapBin` <-
function(object, ...)
{

######################################################################
    ans <- c(list(formula=object$formula1, coefficients=object$results,WALD=object$WALD, AIC=object$aic, empiricalResults=object$empiricalResults, summary.coefs=object$summary.coefs,effect=object$effect))

    class(ans) <- "summary.hapBin"
    return(ans)
}

