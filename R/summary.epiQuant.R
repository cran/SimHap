`summary.epiQuant` <-
function(object, ...)
{
    object.lm <- object$fit.lm
    summary.object.lm <- summary(object.lm)
    ans <- c(summary.object.lm[c("call", "terms", "coefficients", "df", "sigma", "residuals")], list(df.residual=object.lm$df.residual, formula=object$formula, rsquared=summary.object.lm$r.squared, AIC=object$AIC))
    class(ans) <- "summary.epiQuant"
    return(ans)
}
