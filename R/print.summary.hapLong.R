`print.summary.hapLong` <- 
function (x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nLinear mixed-effects model fit by REML\n")
    cat("\nFixed:\n")
    cat(paste(deparse(x$fixed_formula), sep = "\n", collapse = "\n"), sep = "")
    cat("\n\nRandom:\n")
    cat(paste(deparse(x$random_formula), sep = "\n", collapse = "\n"), sep = "")
    
    cat("\n\nHaplotypic effect: ", x$effect, "\n", sep="")
    


        cat("\nCoefficients:\n")
	coefs <- x$summary.coefs
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", has.Pvalue=T, ...)
#####################################	
        
    cat("\n\nAIC: ", x$AIC, "\n", sep = "")

    cat("\n")
    invisible(x)
}
