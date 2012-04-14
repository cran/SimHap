`print.summary.epiClogit` <- 
function (x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
  
    cat("\nCall:\n")
    cat(paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
    resid <- x$residuals
    cat(if (!is.null(x$weights) && diff(range(x$weights)))
        "Weighted ", "Residuals:\n", sep = "")
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- if (length(dim(resid)) == 2)
            structure(apply(t(resid), 1, quantile), dimnames = list(nam,
                dimnames(resid)[[2]]))
        else structure(quantile(resid), names = nam)
        print(rq, digits = digits, ...)

        cat("\nCoefficients:\n")
        coefs <- x$coefficients
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
            na.print = "NA", has.Pvalue=T, ...)

	    
	cat("\nR-squared:\n")
	print(x$rsquared)

    
    cat("\n\nWald Statistic: \n")
    print(x$Wald)

    cat("\n")
    invisible(x)
}
