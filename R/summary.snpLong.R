`summary.snpLong` <- 
function(object, ...)
{
    object.lme <- object$fit.lme
    summary.object.lme <- summary(object.lme)
    resd <- resid(object, type = "pearson")
    if (length(resd) > 5) {
        resd <- quantile(resd)
        names(resd) <- c("Min", "Q1", "Med", "Q3", "Max")
    }

######################################################################
    ans <- c(summary.object.lme[c("call", "terms", "residuals", "fixDF", "sigma", "groups", "dims", "method")], list(coefficients = object$results, random_formula=object$random_formula, fixed_formula=object$fixed_formula, AIC=object$AIC, corStruct=object$corStruct, modelStruct=object.lme$modelStruct))
    
    class(ans) <- "summary.snpLong"
    return(ans)
}
