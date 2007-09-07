`print.summary.hapSurv` <- 
function (x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCox Proportional Hazards Model\n")
    cat("\nCall:\n")
    cat(paste(deparse(x$formula), sep = "\n", collapse = "\n"), sep = "")
    
    cat("\n\nHaplotypic effect: ", x$effect, "\n", sep="")
    
    cat("\nLikelihood Ratio Test: Model without genetic covariates vs model with genetic covariates\n")
    print(x$LRT)
    cat("\n")
    


        cat("\nCoefficients:\n")
	coefs <- x$summary.coefs
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", has.Pvalue=TRUE, ...)
#####################################	
        cat("\nR-squared:\n")
	print(x$rsquared)

    
        cat("\n\nWald Statistic: \n")
        print(x$Wald)

    cat("\n")
    invisible(x)
}
