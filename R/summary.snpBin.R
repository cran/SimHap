`summary.snpBin` <-
function(object, ...)
{
    object.glm <- object$fit.glm
    summary.object.glm <- summary(object.glm)
######################################################################
    ans <- c(summary.object.glm[c("call", "terms","df.residual", "df")], list(residuals=object.glm$residuals,coefficients = object$results, formula=object$formula1, LRT=object$LRT, weights=object.glm$weights, AIC=object$AIC))
    
    class(ans) <- "summary.snpBin"
    return(ans)
}

