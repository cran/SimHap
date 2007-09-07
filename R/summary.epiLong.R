`summary.epiLong` <- 
function(object, ...)
{
    object.lme <- object$fit.lme
    summary.object.lme <- summary(object.lme)
    
    ans <- c(summary.object.lme[c("call", "terms", "fixDF",  "sigma", "groups", "dims", "method", "residuals")], list(coefficients = object$results, fixed_formula=object$fixed_formula, random_formula=object$random_formula, AIC=object$AIC, corStruct=object$corStruct, modelStruct=object.lme$modelStruct))
    
    class(ans) <- "summary.epiLong"
    return(ans)
}

