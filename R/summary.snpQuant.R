`summary.snpQuant` <-
function(object, ...)
{
    object.lm <- object$fit.lm
    summary.object.lm <- summary(object.lm)
######################################################################
    ans <- c(summary.object.lm[c("call", "terms", "na.action", "df", "sigma", "residuals")], list(df.residual=object.lm$df.residual, coefficients = object$results, formula=object$formula1, LRT=object$LRT, weights=object.lm$weights, AIC=object$AIC, rsquared=object$rsquared, predicted.values=object$predicted.values))
    
    class(ans) <- "summary.snpQuant"
    return(ans)
}

