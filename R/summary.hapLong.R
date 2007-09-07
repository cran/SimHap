`summary.hapLong` <- 
function(object, ...)
{
    
######################################################################
    ans <- c(list(fixed_formula=object$fixed_formula, random_formula=object$random_formula, coefficients=object$results, empiricalResults=object$empiricalResults, results=object$results, summary.coefs=object$summary.coefs, AIC=object$AIC,corStruct=object$corStruct, effect=object$effect))
    
    class(ans) <- "summary.hapLong"
    return(ans)
}
