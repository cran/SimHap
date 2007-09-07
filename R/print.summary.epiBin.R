`print.summary.epiBin` <-
function (x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
  
    cat("\nCall:\n")
    cat(paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
    resid <- x$residuals
    df <- x$df
    cat(if (!is.null(x$weights) && diff(range(x$weights)))
        "Weighted ", "Residuals:\n", sep = "")
    if (x$df[2] > 5) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- if (length(dim(resid)) == 2)
            structure(apply(t(resid), 1, quantile), dimnames = list(nam,
                dimnames(resid)[[2]]))
        else structure(quantile(resid), names = nam)
        print(rq, digits = digits, ...)
    }
    else if (x$df[2] > 0) {
        print(resid, digits = digits, ...)
    }
    else {
        cat("ALL", df[1], "residuals are 0: no residual degrees of freedom!\n")
    }
    
        if (!is.null(df <- x$df) && (nsingular <- df[3] - df[1]))
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients

        printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
            na.print = "NA", has.Pvalue=TRUE, ...)
    
    cat("\n\nAIC: ", format(x$AIC, digits = max(4, digits + 1)), "\n", sep = "")

    cat("\n")
    invisible(x)
}

