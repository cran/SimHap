`print.summary.hapQuant` <-
function (x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
  
    cat("\nCall:\n")
    cat(paste(deparse(x$formula), sep = "\n", collapse = "\n"), sep = "")
    cat("\nHaplotypic effect: ", x$effect, "\n\n", sep="")
    
    cat("\nLikelihood Ratio Test: Model without genetic covariates vs model with genetic covariates\n")
    print(x$LRT)
    cat("\n")
    

        cat("\nCoefficients:\n")
        coefs <- x$summary.coefs
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", has.Pvalue=TRUE, ...)
#####################################
        cat("\nEstimated marginal means of dependent variable broken down by haplotype levels, evaluated at mean values of the model predictors:\n")
print(x$predicted) 

        cat("\nAdjusted R-squared:\n")
print(x$rsquared)

    cat("\n\nAIC: ", x$AIC, "\n", sep = "")

    cat("\n")
    invisible(x)
}

