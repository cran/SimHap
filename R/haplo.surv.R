`haplo.surv` <-
function(formula1, formula2, pheno, haplo, sim, effect="add", sub=NULL) {

  hapFreqs <- haplo$hapObject$final.freq
  haplo <- haplo$hapData
  library(stats)
  call <- match.call()
  library(survival)

  formula1 <- formula(formula1)
  formula2 <- formula(formula2)

  formula1_nofactors <- formula1
  formula1_terms <- attr(terms(formula1_nofactors), "term.labels")

  if(any(regexpr(":", formula1_terms)!=-1)){
        formula1_terms <- formula1_terms[-which(regexpr(":", formula1_terms)!=-1)]
  }

  if(any(regexpr("factor", formula1_terms)==1)) {
        formula1_terms[which(regexpr("factor", formula1_terms)==1)] <- substr(formula1_terms[which(regexpr("factor", formula1_terms)==1)],8,nchar(formula1_terms[which(regexpr("factor", formula1_terms)==1)])-1)
  }
  #else formula1_terms <- attr(terms(formula1_nofactors), "term.labels")

  freq.estnums <- freqTest(terms=formula1_terms, freqs=hapFreqs, n=length(unique(haplo[,1])), effect=effect)

  # first column retains number of non-zero weights for individual
  # second column holds the current iteration index for the next weight change
  num_weights <- matrix(0, nrow=nrow(pheno), ncol=2)
  num_indivs <- length(unique(as.numeric(haplo[,1])))
  #num_indivs <- as.numeric(haplo[nrow(haplo),1])

  # these are the two distributed occurrence matrices representing both corresponding haplotypes
  # the dimensions are #iterations by #individuals
  hap1s_result <- matrix(0, nrow=sim, ncol=num_indivs)
  hap2s_result <- matrix(0, nrow=sim, ncol=num_indivs)

  # *****************************
  print("* Finding highest individual frequency ...")
  # this section calculates the highest frequency of an indiv

  lastID <- haplo[1,1]
  count <- 1
  biggest <- 1

  for(i in 2:nrow(haplo)) {
    tmpID <- haplo[i,1]

    if(lastID==tmpID) {

      # only increment count if the weight is not too small in the context
      # of the number of iterations
        count <- count+1
    }
    else {
      if(count>biggest)
        biggest <- count

      count <- 1
    }

    lastID <- tmpID
  }

  # at this point 'biggest' is the highest frequency of an individual
  # which translates to the highest number of different haplotype combinations
  print("  Done")
  # *****************************

  indiv_weights <- matrix(0, nrow=num_indivs, ncol=biggest)
  indiv_hap1s <- matrix(0, nrow=num_indivs, ncol=biggest)
  indiv_hap2s <- matrix(0, nrow=num_indivs, ncol=biggest)


  lastID <- haplo[1,1]
  indiv_weights[1,1] <- haplo[1,4]
  indiv_hap1s[1,1] <- haplo[1,2]
  indiv_hap2s[1,1] <- haplo[1,3]
  count <- 1
  indiv <- 1

  # ****************************
  print("* Populating individual haplotypes and posterior probabilities ...")
  # This section makes a count of the number of occurences of each individual

  for(i in 2:nrow(haplo)) {
    tmpID <- haplo[i,1]
    this_weight <- as.numeric(haplo[i,ncol(haplo)])

    # one attempt at rounding too small weight*sim up to '1'
    if((sim*this_weight) < 1)
      this_weight <- 1 / sim

    # the next element is the same ID as the last
    if(lastID==tmpID) {

      # only increment count if the weight is not too small in the context
      # of the number of iterations
      if((sim*this_weight) >= 1)
        count <- count+1
    }
    # this element's ID differs from the last: store and restart count
    else {
      num_weights[indiv,1] <- count
      indiv <- indiv + 1
      count <- 1
    }

    if((sim*this_weight) >= 1) {
      indiv_weights[indiv,count] <- this_weight
      indiv_hap1s[indiv,count] <- haplo[i,2]
      indiv_hap2s[indiv,count] <- haplo[i,3]
    }

    lastID <- tmpID
  }

  # take care of the last element
  num_weights[indiv,1] <- count

  print("  Done")
  # ****************************

  # ****************************
  print("* Distributing individual occurrences across the simulations by posterior probability ...")
  # main loop across all iterations
  if(length(sub)==0) {
      fit2.coxph <- coxph(formula=formula2, data=pheno, na.action=na.omit)
      }

     else {
      fit2.coxph <- eval(substitute(coxph(formula=formula2, data=pheno, na.action=na.omit, subset=subset), list(subset=sub)))
      }

  fit2.rsquared <- summary(fit2.coxph)$rsq[1]
  lnLsmall <- fit2.coxph$loglik[2]
  fit2.df <- summary(fit2.coxph)$logtest[2]

  for(i in 1:sim) {

    # determine weight for each individual and populate hapXs vectors
    for(j in 1:num_indivs) {

      # save processing time ... ;)
      current_numw <- num_weights[j,1]

      # invalid data case, weights for an indiv do not add up to one
      if(current_numw == 0) {
        print("Error. Weights do not sum to 1, indiv ID:")
        print(j)
        print(i)
        print("NOTE: Intermediate result returned")
        return(weights.result)
      }

      weight <- as.numeric(indiv_weights[j,current_numw])
      hap1s_result[i,j] <- indiv_hap1s[j,current_numw]
      hap2s_result[i,j] <- indiv_hap2s[j,current_numw]

      if(i>=(num_weights[j,2] + (sim*weight))) {
        num_weights[j,1] <- current_numw - 1
        num_weights[j,2] <- num_weights[j,2] + (sim*weight)
      }
    }

    # report on progress
    if(i==round(sim*0.01))
      print("    1%")
    if(i==round(sim*0.05))
      print("    5%")
    if(i==round(sim*0.25))
      print("    25%")
    if(i==round(sim*0.5))
      print("    50%")
    if(i==round(sim*0.75))
      print("    75%")
    if(i==round(sim*0.90))
      print("    90%")

  }

  print("  Done")
  # ****************************

  # ****************************
  print("* Generating a random pattern of individuals for each simulation ...")
  # This section generates the random choice of individual

  # create a vector with a linear progression 0 -> sim
  choice <- 1:sim
  sim_choice <- matrix(0, nrow=num_indivs, ncol=sim)

  # generate a random choice for each individual
  for(i in 1:num_indivs) {
    sim_choice[i,] <- sample(choice, sim, replace=F)
  }

  print("  Done")
  # ****************************


  # ****************************
  print("* Constructing dataframes and performing cox proportional hazards model for each simulation ...")
  # This section constructs the dataframe for each iteration in preparation for coxph

  haplo_table <- table(c(haplo[,2],haplo[,3]))
  num_haplos <- dim(haplo_table)
  names_haplos <- names(haplo_table)

  # prepare the reusable dataframe container ...
  dataframe_extra <- matrix(0, nrow=num_indivs, ncol=num_haplos)

  # perform loop through all iterations, constructing the dataframe and applying coxph model to each one
#------------------------------------

    coxph.out <- function(x, digits=max(options()$digits - 4, 3), ...)
    {

      if (!is.null(x$fail)) {
        cat(" Coxph failed.", x$fail, "\n")
        return()
        }

      savedig <- options(digits = digits)
      on.exit(options(savedig))

      coef <- x$coef
      se <- sqrt(diag(x$var))
      if(is.null(coef) | is.null(se))
        stop("Input is not valid")

      if (is.null(x$naive.var)) {
        tmp <- cbind(coef, exp(coef), se, coef/se,
        signif(1 - pchisq((coef/ se)^2, 1), digits -1))
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "z", "p"))
        }

      else {
        nse <- sqrt(diag(x$naive.var))
        tmp <- cbind(coef, exp(coef), nse, se, coef/se,
        signif(1 - pchisq((coef/se)^2, 1), digits -1))
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "robust se", "z", "p"))
        }
      tmp <- as.data.frame(tmp)
      return(tmp)
    }


  coef.dat <- NULL
  p.dat <- NULL
  stderror.dat <- NULL
  OR.dat <- NULL
  ORupper.dat <- NULL
  ORlower.dat <- NULL
  out <- NULL
  output <- NULL
  fit1.rsquared.dat <- NULL
  max.rsquared.dat <- NULL
  stderror.dat <- NULL

  anovp.dat <- NULL
  anovdf.dat <- NULL
  anovloglik.dat <- NULL
  anov.out <- NULL
  anov.out1 <- NULL
  anov.out2 <- NULL
  lnLbig.dat <- NULL
  lr.dat <- NULL
  lrt.dat <- NULL

  # the dynamic point at which the algorithm reports progress
  five_percent <- round(sim*0.05)
  report <- five_percent

  # main loop
  for(i in 1:sim) {
    dataframe_extra[,] <- 0

    # Dominant model
    if(effect=="dom") {
      for(j in 1:num_indivs) {
        simul <- sim_choice[j,i]
        hap1_str <- hap1s_result[simul, j]
        hap2_str <- hap2s_result[simul, j]

        colnumber <- match(hap1_str, names_haplos)
        dataframe_extra[j,colnumber] <- 1

        if(hap1_str != hap2_str) {
          colnumber <- match(hap2_str, names_haplos)
          dataframe_extra[j,colnumber] <- 1
        }


      }
    }
    # Recessive model
    if(effect=="rec") {
      for(j in 1:num_indivs) {
        simul <- sim_choice[j,i]
        hap1_str <- hap1s_result[simul, j]

        if(hap1_str == hap2s_result[simul, j]) {
          colnumber <- match(hap1_str, names_haplos)
          dataframe_extra[j,colnumber] <- 1
        }


      }
    }

    # Additive model
    if(effect=="add") {
      for(j in 1:num_indivs) {
        simul <- sim_choice[j,i]
        colnumber <- match(hap1s_result[simul, j], names_haplos)
        dataframe_extra[j,colnumber] <- 1
        colnumber <- match(hap2s_result[simul, j], names_haplos)
        dataframe_extra[j,colnumber] <- dataframe_extra[j,colnumber] + 1
      }
    }

    # Expand dataframe_extra to match dimension of data for sake of R argument
    dataframe_extra_expand <- matrix(0, nrow=nrow(pheno), ncol=ncol(dataframe_extra))


    de_row <- 1

    dataframe_extra_expand[1,] <- dataframe_extra[de_row,]
    last_phen_ID <- pheno[1,1]


    phen_ID <- pheno[,1]
    for(h in 2:nrow(pheno)) {
      current_phen_ID <- phen_ID[h]

      if (current_phen_ID != last_phen_ID) {
        de_row <- de_row + 1
      }



      dataframe_extra_expand[h,] <- dataframe_extra[de_row,]

      last_phen_ID <- current_phen_ID
    }

    # concatenate the extra haplotype columns
    dataframe_extra_expand <- as.data.frame(dataframe_extra_expand)
    for(j in 1:ncol(dataframe_extra_expand)){
        colnames(dataframe_extra_expand)[j] <- names_haplos[j]
        }

    dataframe <- as.data.frame(cbind(pheno, dataframe_extra_expand))
    dataframe <- dataframe[complete.cases(dataframe[formula1_terms]),]

    # perform the coxph with the current dataframe
    # coxph

     if(length(sub)==0) {
      fit.coxph <- coxph(formula=formula1, data=dataframe, na.action=na.omit)
      }

     else {
      fit.coxph <- eval(substitute(coxph(formula=formula1, data=dataframe, na.action=na.omit, subset=subset), list(subset=sub)))
      }

    fit1.rsquared <- summary(fit.coxph)$rsq[1]
    fit1.rsquared.dat <- rbind(fit1.rsquared.dat, fit1.rsquared)
    max.rsquared <- summary(fit.coxph)$rsq[2]
    max.rsquared.dat <- rbind(max.rsquared.dat, max.rsquared)

    lnLbig <- fit.coxph$loglik[2]
    fit1.df <- summary(fit.coxph)$logtest[2]
    lnLbig.dat <- rbind(lnLbig.dat, lnLbig)
    fit2.df <- summary(fit2.coxph)$logtest[2]

    lr <- -2*(lnLsmall-lnLbig)
    lr.dat <- rbind(lr.dat, lr)
    lr.df <- fit1.df-fit2.df
    lrt <- pchisq(lr,df=lr.df)
    lrt.dat <- rbind(lrt.dat, lrt)

    anov <- as.data.frame(anova(fit2.coxph, fit.coxph, test="Chisq"))

    wald <- formatC(summary(fit.coxph)$waldtest)
    wald <- as.data.frame(cbind(wald[1], wald[2], wald[3]))
    row.names(wald) <- ""
    names(wald) <- c("Wald test", "df", "P-value")

    out.coxph <- coxph.out(fit.coxph)
    out.coxph$ORlower <- exp(out.coxph$coef - 1.96 * out.coxph$se)
    out.coxph$ORupper <- exp(out.coxph$coef + 1.96 * out.coxph$se)

    # extract some elements from the coxph fitted model
    # add this row to coef.dat

    anovp.row <- anov$"P(>|Chi|)"[2]
    anovp.dat <- rbind(anovp.dat, anovp.row)
    anovdf.row <- anov$Df
    anovdf.dat <- rbind(anovdf.dat, anovdf.row)
    anovloglik.row <- anov$loglik
    anovloglik.dat <- rbind(anovloglik.dat, anovloglik.row)
    coef.row <- out.coxph$coef
    coef.dat <- rbind(coef.dat, coef.row)
    pval.row <- out.coxph$p
    p.dat <- rbind(p.dat, pval.row)
    stderror.row <- out.coxph$"se(coef)"
    stderror.dat <- rbind(stderror.dat, stderror.row)

    # add a row to stderror.dat

    OR.row <- out.coxph$"exp(coef)"
    OR.dat <- rbind(OR.dat, OR.row)

    ORupper.row <- out.coxph$ORupper
    ORupper.dat <- rbind(ORupper.dat, ORupper.row)

    ORlower.row <- out.coxph$ORlower
    ORlower.dat <- rbind(ORlower.dat, ORlower.row)

    # report on progress
    if(i==report) {
      percentage <- report * 100 / sim
      print(paste(percentage, "%"))
      report <- report + five_percent
    }
  }

  # nullify the row names

  row.names(anovdf.dat) <- NULL
  row.names(anovp.dat) <- NULL
  row.names(anovloglik.dat) <- NULL
  row.names(coef.dat) <- NULL
  row.names(stderror.dat) <- NULL
  row.names(p.dat) <- NULL
  row.names(OR.dat) <- NULL
  row.names(ORupper.dat) <- NULL
  row.names(ORlower.dat) <- NULL
  row.names(lnLbig.dat) <- NULL
  row.names(fit1.rsquared.dat) <- NULL
  row.names(max.rsquared.dat) <- NULL

  coef.dat <- as.data.frame(coef.dat)
  stderror.dat <- as.data.frame(stderror.dat)
  p.dat <- as.data.frame(p.dat)
  OR.dat <- as.data.frame(OR.dat)
  ORupper.dat <- as.data.frame(ORupper.dat)
  ORlower.dat <- as.data.frame(ORlower.dat)

  anovdf.dat <- as.data.frame(anovdf.dat)
  anovp.dat <- as.data.frame(anovp.dat)
  anovloglik.dat <- as.data.frame(anovloglik.dat)
  lnLbig.dat <- as.data.frame(lnLbig.dat)
  lr.dat <- as.data.frame(lr.dat)
  lrt.dat <- as.data.frame(lrt.dat)

  allResults <- list(HR=OR.dat, HRupper=ORupper.dat, HRlower=ORlower.dat, P.Value=p.dat)
  names(allResults$HR) <- row.names(out.coxph)
  names(allResults$HRupper) <- row.names(out.coxph)
  names(allResults$HRlower) <- row.names(out.coxph)
  names(allResults$P.Value) <- row.names(out.coxph)

  sum.of.squares <- NULL
  for(i in 1:ncol(stderror.dat)){
       sum.of.squares <- cbind(sum.of.squares,sum(stderror.dat[,i]^2))
  }
  sum.of.squares <- as.data.frame(sum.of.squares)
  names(sum.of.squares) <- names(stderror.dat)
  se1 <- sqrt(sum.of.squares/nrow(stderror.dat))
  se2 <- sapply(coef.dat, sd)
  se.adj <- sqrt(se1^2 + se2^2)

  output$OR <- formatC(colMeans(OR.dat))
  output$ORupper <- formatC(colMeans(ORupper.dat))
  output$ORlower <- formatC(colMeans(ORlower.dat))
  output$pval <- formatC(colMeans(p.dat))

  summary.coefs <- data.frame(cbind(colMeans(allResults$HR), exp(colMeans(coef.dat)-(1.96*as.numeric(as.vector(se.adj)))), exp(colMeans(coef.dat)+(1.96*as.numeric(as.vector(se.adj)))), colMeans(allResults$P.Value)))
  names(summary.coefs) <- c("Hazard.Ratio", "HR.lower", "HR.upper", "P.Value")

  LRT.out1 <- cbind(round(colMeans(lnLbig.dat),digits=4), fit1.df, round(colMeans(lr.dat),digits=4), signif((1-colMeans(lrt.dat)), digits=4))
  LRT.out2 <- cbind(round(lnLsmall,digits=4), fit2.df, "", "")
  row.names(LRT.out1) <- c("Full model")
  row.names(LRT.out2) <- c("Non-genetic")
  LRT <- rbind(LRT.out1, LRT.out2)
  LRT <- as.data.frame(LRT)
  names(LRT) <- c("logLik", "df", "LR", "P.Value")

  anov.out1 <- cbind(round(colMeans(anovloglik.dat[1]), digits=3), "", "")
  row.names(anov.out1) <- c("1")
  anov.out2 <- cbind(round(colMeans(anovloglik.dat[2]), digits=3), colMeans(anovdf.dat[2]), signif(colMeans(anovp.dat), digits=3))
  row.names(anov.out2) <- c("2")
  likelihood.out <- paste("'log Lik'", round(colMeans(lnLbig.dat), digits=3), paste("(df=", fit1.df, ")", sep=""))

  rsquared.out <- rbind(round(fit2.rsquared,digits=3),round(colMeans(fit1.rsquared.dat),digits=3), round(colMeans(max.rsquared.dat), digits=3))
  rsquared.out <- as.data.frame(rsquared.out)
  names(rsquared.out) <- "R-Squared"
  row.names(rsquared.out) <- c("Without Haplotypes", "Including Haplotypes", "Max R-Squared")

  print("  Done")
  # ****************************

  # Arrange the output data

  for(i in 1:ncol(coef.dat)){
    output$OR.CI[i] <- paste("(",formatC(quantile(OR.dat[,i], probs=c(0.025), na.rm=T)),",",formatC(quantile(OR.dat[,i], probs=c(0.975), na.rm=T)),")", sep="")
    output$ORlower.CI[i] <- paste("(",formatC(quantile(ORlower.dat[,i], probs=c(0.025), na.rm=T)),",",formatC(quantile(ORlower.dat[,i], probs=c(0.975), na.rm=T)),")", sep="")
    output$ORupper.CI[i] <- paste("(",formatC(quantile(ORupper.dat[,i], probs=c(0.025), na.rm=T)),",",formatC(quantile(ORupper.dat[,i], probs=c(0.975), na.rm=T)),")", sep="")
    output$pval.CI[i] <- paste("(",formatC(quantile(p.dat[,i], probs=c(0.025), na.rm=T)),",",formatC(quantile(p.dat[,i], probs=c(0.975), na.rm=T)),")", sep="")
#    output$se.CI[i] <- paste("(",formatC(quantile(stderror.dat[,i], probs=c(0.025), na.rm=T)),",",formatC(quantile(stderror.dat[,i], probs=c(0.975), na.rm=T)),")", sep="")
  }


  out <- data.frame(cbind(output$OR, output$OR.CI, output$ORlower, output$ORlower.CI, output$ORupper, output$ORupper.CI, output$pval, output$pval.CI))

  names(out) <- c("HR", "HR.quantiles", "HRlower.95", "HRlower.quantiles", "HRupper.95", "HRupper.quantiles", "P.Val", "P.Val.quantiles")
  row.names(out) <- row.names(out.coxph)
  anov.out <- rbind(anov.out1, anov.out2)
  anov.out <- as.data.frame(anov.out)

  names(anov.out) <- c("Residual DF", "DF", "P-Value")

  if(effect=="add") Effect <- ("ADDITIVE")
  if(effect=="dom") Effect <- ("DOMINANT")
  if(effect=="rec") Effect <- ("RECESSIVE")

  out.list <- list(formula1=formula1, formula2=formula2, results=out,empiricalResults=allResults, summary.coefs=summary.coefs, logLik=likelihood.out,LRT=LRT, ANOVA=anov.out,Wald=wald, rsquared=rsquared.out, effect=Effect)
  class(out.list) <- "hapSurv"
  return(out.list)

  return(out.list)
}

