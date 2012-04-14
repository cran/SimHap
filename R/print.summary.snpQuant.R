`print.summary.snpQuant` <-
function (x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
  
    x.lm <- x$fit.lm
    cat("\nCall:\n")
    cat(paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
    cat("Likelihood Ratio Test: Model without genetic covariates vs model with genetic covariates:\n")
    print(x$LRT)
    cat("\n")
    
    #cat("Residuals: \n")
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
    
        if (!is.null(df <- x.lm$df) && (nsingular <- df[3] - df[1]))
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients

        printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
            na.print = "NA", has.Pvalue=TRUE, ...)
    
cat("\nEstimated marginal means of dependent variable broken down by genotype, evaluated at mean values of the model predictors:\n")
print(x$predicted.values)
#####################################    
        cat("\nResidual standard error:", format(signif(x$sigma,
        digits)), "on", x$df[2], "degrees of freedom\n")
    if (nchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")
    #if (!is.null(x$fstatistic)) {
cat("\nAdjusted R-squared:\n")
print(x$rsquared)
    #}
######################################
   
    #cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
     #   format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null",
     #       "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance",
     #       "deviance")]), digits = max(5, digits + 1)), " on",
     #       format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
     #       1, paste, collapse = " "), sep = "")
    
    cat("\n\nAIC: ", format(x$AIC, digits = max(4, digits + 1)), "\n", sep = "")

    cat("\n")
    invisible(x)
}

