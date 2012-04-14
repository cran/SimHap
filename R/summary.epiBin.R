`summary.epiBin` <-
function(object, ...)
{
    object.glm <- object$fit.glm
    summary.object.glm <- summary(object.glm)
######################################################################
    ans <- c(summary.object.glm[c("call", "terms","df.residual", "df")], list(coefficients=object$results,residuals=object.glm$residuals,formula=object$formula, weights=object.glm$weights, AIC=object$AIC))
    
    class(ans) <- "summary.epiBin"
    return(ans)
}
