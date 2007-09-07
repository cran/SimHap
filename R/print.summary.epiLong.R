`print.summary.epiLong` <- 
function (x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
    dd <- x$dims
    
    cat("Linear mixed-effects model fit by ")
    cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))

    print(summary(x$modelStruct), sigma = x$sigma, reEstimates = x$coef$random)
    cat("\nFixed effects: ")
    print(x$fixed_formula)
    
    cat("\nCoefficients:\n")
    coefs <- x$coefficients

    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", has.Pvalue=T, ...)    

    cat("\nStandardized Within-Group Residuals:\n")
    print(x$residuals)
    cat("\nNumber of Observations:", x$dims[["N"]])
    cat("\nNumber of Groups: ")
    Ngrps <- dd$ngrps[1:dd$Q]
    if ((lNgrps <- length(Ngrps)) == 1) {
        cat(Ngrps, "\n")
    }
    else {
        sNgrps <- 1:lNgrps
        aux <- rep(names(Ngrps), sNgrps)
        aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps,
            lNgrps))[!lower.tri(diag(lNgrps))])
        names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
        cat("\n")
        print(rev(Ngrps))
    }
  
    cat("\n\nAIC: ", format(x$AIC, digits = max(4, digits + 1)), "\n", sep = "")

    cat("\n")
    invisible(x)
}
