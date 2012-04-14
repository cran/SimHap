`haplo.bin` <-
function(formula1, formula2, pheno, haplo, sim, effect="add", sub=NULL, adjust=FALSE) {
  library(stats)
  call <- match.call()

  hapFreqs <- haplo$hapObject$final.freq
  haplo <- haplo$hapData

  if(!identical(as.character(unique(pheno$ID)), as.character(unique(haplo$ID)))) stop("Phenotype data and Haplotype data are not in the same order.")

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
  num_indivs <- nrow(pheno)

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
      #if((sim*(as.numeric(haplo[i,ncol(haplo)]))) > 1)
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
  indiv_hap1s[1,1] <- haplo[1,2]
  indiv_weights[1,1] <- haplo[1,4]
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
  fit2.glm <- eval(substitute(glm(formula2, data=pheno, family=binomial, subset=subset), list(subset=sub)))

  #log likelihood for smaller, nested model (without haplotypes)
  lnLsmall <- logLik(fit2.glm)

  for(i in 1:sim) {

    # determine weight for each individual and populate hapXs vectors
    for(j in 1:num_indivs) {

      # save processing time ... ;)
      current_numw <- num_weights[j,1]

      # invalid pheno case, weights for an indiv do not add up to one
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
    sim_choice[i,] <- sample(choice, sim, replace=FALSE)
  }

  print("  Done")
  # ****************************


  # ****************************
  print("* Constructing dataframes and performing generalised linear model for each simulation ...")
  # This section constructs the dataframe for each iteration in preparation for glm

  haplo_table <- table(c(haplo[,2],haplo[,3]))
  num_haplos <- dim(haplo_table)
  names_haplos <- names(haplo_table)


  # prepare the reusable dataframe container ...
  dataframe_extra <- matrix(0, nrow=num_indivs, ncol=num_haplos)
  dataframe_extra <- as.data.frame(dataframe_extra)
  for(i in 1:ncol(dataframe_extra)){
    colnames(dataframe_extra)[i] <- paste(names_haplos[i])}

  # perform loop through all iterations, constructing the dataframe and applying glm to each one
#------------------------------------

  coef.dat <- NULL
  p.dat <- NULL
  stderror.dat <- NULL
  out <- NULL
  output <- NULL
  pvals <- NULL
  stderrors <- NULL

  anovp.dat <- NULL
  anovdf.dat <- NULL
  anovresdf.dat <- NULL

  anovfullp.dat <- NULL
  anovfulldf.dat <- NULL
  anovfullresdf.dat <- NULL

  anov.out <- NULL
  anov.out1 <- NULL
  anov.out2 <- NULL
  anovfull.out <- NULL

  aic.dat <- NULL
  lr.dat <- NULL
  lrt.dat <- NULL
  lnLbig.dat <- NULL
  vcov.list <- list(NULL)
  beta.list <- list(NULL)

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

    # concatenate the extra haplotype columns
    dataframe <- as.data.frame(cbind(pheno, dataframe_extra))
    #change dataframe (if necessary) to include only indivs with complete data for all terms in formula1
    dataframe <- dataframe[complete.cases(dataframe[formula1_terms]),]

    # perform the glm with the current dataframe
    # glm
    fit1.glm <- eval(substitute(glm(formula1, data=dataframe, family=binomial, subset=subset), list(subset=sub)))

    fit.glm <- as.data.frame(summary(fit1.glm)$coefficients)
    anov <- as.data.frame(anova(fit2.glm, fit1.glm, test="Chisq"))
    anovfull <- as.data.frame(anova(fit1.glm, test="Chisq"))

    #extract log-likelihood of model with haplotypes
    lnLbig <- logLik(fit1.glm)
    lnLbig.dat <- rbind(lnLbig.dat, lnLbig)
    lr <- -2*(lnLsmall[1]-lnLbig[1])
    lr.dat <- rbind(lr.dat, lr)
    lr.df <- attr(lnLbig, "df")-attr(lnLsmall,"df")
    lrt <- pchisq(lr,df=lr.df)
    lrt.dat <- rbind(lrt.dat, lrt)

    # extract variance-covariance matrix
    vcov.list[[i]] <- vcov(fit1.glm)
    beta.list[[i]] <- fit.glm$Estimate

    aic <- AIC(fit1.glm)
    aic.dat <- rbind(aic.dat, aic)

    # add this row to anovfull.dat
    anovfullp.row <- anovfull$"Pr(>Chi)"
    anovfullp.dat <- rbind(anovfullp.dat, anovfullp.row)
    anovfulldf.row <- anovfull$Df
    anovfulldf.dat <- rbind(anovfulldf.dat, anovfulldf.row)
    anovfullresdf.row <- anovfull$"Resid. Df"
    anovfullresdf.dat <- rbind(anovfullresdf.dat, anovfullresdf.row)

    # add this row to anovp.dat

    anovp.row <- anov$"Pr(>Chi)"[2]
    anovp.dat <- rbind(anovp.dat, anovp.row)
    anovdf.row <- anov$Df
    anovdf.dat <- rbind(anovdf.dat, anovdf.row)
    anovresdf.row <- anov$"Resid. Df"
    anovresdf.dat <- rbind(anovresdf.dat, anovresdf.row)

    # add this row to coef.dat
    coef.row <- fit.glm$Estimate
    coef.dat <- rbind(coef.dat, coef.row)
    stderror.row <- fit.glm$"Std. Error"
    stderror.dat <- rbind(stderror.dat, stderror.row)


    # extract some elements from the glm summary method and add row to p.dat
    pvals <- t(fit.glm[ncol(fit.glm)])
    p.dat <- rbind(p.dat, pvals)

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
  row.names(anovresdf.dat) <- NULL

  row.names(anovfulldf.dat) <- NULL
  row.names(anovfullp.dat) <- NULL
  row.names(anovfullresdf.dat) <- NULL

  row.names(p.dat) <- NULL
  row.names(stderror.dat) <- NULL
  row.names(coef.dat) <- NULL

  row.names(aic.dat) <- NULL
  row.names(lr.dat) <- NULL
  row.names(lrt.dat) <- NULL
  row.names(lnLbig.dat) <- NULL

  anovdf.dat <- as.data.frame(anovdf.dat)
  anovp.dat <- as.data.frame(anovp.dat)
  anovresdf.dat <- as.data.frame(anovresdf.dat)

  anovfulldf.dat <- as.data.frame(anovfulldf.dat)
  anovfullp.dat <- as.data.frame(anovfullp.dat)
  anovfullresdf.dat <- as.data.frame(anovfullresdf.dat)

  aic.dat <- as.data.frame(aic.dat)
  lr.dat <- as.data.frame(lr.dat)
  lrt.dat <- as.data.frame(lrt.dat)
  lnLbig.dat <- as.data.frame(lnLbig.dat)

  p.dat <- as.data.frame(p.dat)
  coef.dat <- as.data.frame(coef.dat)
  stderror.dat <- as.data.frame(stderror.dat)

  names(aic.dat) <- c("AIC")

  allResults <- list(Coef=coef.dat, Std.Error=stderror.dat, P.Value=p.dat)
  names(allResults$Coef) <- row.names(fit.glm)
  names(allResults$Std.Error) <- row.names(fit.glm)
  names(allResults$P.Value) <- row.names(fit.glm)

#  allResults <- list(OR=OR.dat, OR.lower.95CI=ORlower.dat, OR.upper.95CI=ORupper.dat, P.Value=p.dat)
#  names(allResults$OR) <- row.names(fit.glm)
#  names(allResults$OR.lower.95CI) <- row.names(fit.glm)
#  names(allResults$OR.upper.95CI) <- row.names(fit.glm)
#  names(allResults$P.Value) <- row.names(fit.glm)

#  sum.of.squares <- NULL
#  for(i in 1:ncol(stderror.dat)){
#       sum.of.squares <- cbind(sum.of.squares,sum(stderror.dat[,i]^2))
#  }
#  sum.of.squares <- as.data.frame(sum.of.squares)
#  names(sum.of.squares) <- names(stderror.dat)
#  se1 <- sqrt(sum.of.squares/nrow(stderror.dat))
#  se2 <- sd(coef.dat)
#  se.adj <- sqrt(se1^2 + se2^2)

  # Combine inferences across the imputed datasets
  out.mi <- UVI(coef.dat, stderror.dat^2,n=num_indivs, ADJ=adjust)

  ind.haploeffect <- which(!is.element(names(fit1.glm$coefficients), names(fit2.glm$coefficients)))
  p.full <- length(fit1.glm$coefficients)
  L.contrast <- NULL
  for(j in 1:length(ind.haploeffect)){
    L.contrast <- rbind(L.contrast, c(rep(0, ind.haploeffect[j]-1),1,rep(0, p.full-ind.haploeffect[j])) )
  }
  out.mi.haps <- MVI(beta.list, vcov.list, L=L.contrast)
  out.mi.haps <- out.mi.haps
  out.coef <- as.numeric(formatC(out.mi$coefficients))
  out.pval <- as.numeric(formatC(out.mi$p.value))
  out.se <- as.numeric(formatC(out.mi$se))
  #if(!is.null(predicted.dat)) predicted.vals <- formatC(mean(predicted.dat))

  summary.coefs <- data.frame(cbind(out.coef, out.se, out.pval), row.names=row.names(fit.glm))
  names(summary.coefs) <- c("Coefficient", "Std.error", "P.Value")

  WALD.out <- cbind(round(out.mi.haps[4]), round(out.mi.haps[5],2), round(out.mi.haps[1], digits=4), round(out.mi.haps[3], digits=4))
  WALD.out <- as.data.frame(WALD.out)
  names(WALD.out) <- c("Num DF","Den DF","F.Stat", "P.Value")
  row.names(WALD.out) <- ""


  anovfull.out <- cbind(colMeans(anovfullresdf.dat), colMeans(anovfulldf.dat), formatC(colMeans(anovfullp.dat)))
  row.names(anovfull.out) <- row.names(anovfull)
  anovfull.out[1,2] <- ""
  anovfull.out[1,3] <- ""
  anovfull.out <- as.data.frame(anovfull.out)
  names(anovfull.out) <- c("Residual DF", "DF", "P-Value")
  likelihood.out <- paste("'log Lik'", round(colMeans(lnLbig.dat), digits=3), paste("(df=", attr(lnLbig, "df"), ")", sep=""))

  anov.out1 <- cbind(colMeans(anovresdf.dat[1]), "", "")
  row.names(anov.out1) <- c("1")
  anov.out2 <- cbind(colMeans(anovresdf.dat[2]), colMeans(anovdf.dat[2]), signif(colMeans(anovp.dat), digits=3))
  row.names(anov.out2) <- c("2")



  print("  Done")
  # ****************************

  # Arrange the output data

  for(i in 1:ncol(coef.dat)){

    out$coef.CI[i] <- paste("(",formatC(quantile(coef.dat[,i], probs=c(0.025), na.rm=TRUE)),",",formatC(quantile(coef.dat[,i], probs=c(0.975), na.rm=TRUE)),")", sep="")
    out$pval.CI[i] <- paste("(",formatC(quantile(p.dat[,i], probs=c(0.025), na.rm=TRUE)),",",formatC(quantile(p.dat[,i], probs=c(0.975), na.rm=TRUE)),")", sep="")
    out$se.CI[i] <- paste("(",formatC(quantile(stderror.dat[,i], probs=c(0.025), na.rm=TRUE)),",",formatC(quantile(stderror.dat[,i], probs=c(0.975), na.rm=TRUE)),")", sep="")

  }
  out <- data.frame(cbind(out.coef, out$coef.CI, out.se, out$se.CI, out.pval, out$pval.CI))
  names(out) <- c("Coef", "Coef.quantiles", "Std.Error", "Std.Error.quantiles", "P.Val", "P.Val.quantiles")
  row.names(out) <- row.names(fit.glm)


  anov.out <- rbind(anov.out1, anov.out2)
  anov.out <- as.data.frame(anov.out)

  names(anov.out) <- c("Residual DF", "DF", "P.Value")

  if(effect=="add") Effect <- ("ADDITIVE")
  if(effect=="dom") Effect <- ("DOMINANT")
  if(effect=="rec") Effect <- ("RECESSIVE")

  "%w/o%" <- function(x,y) x[!x %in% y]
  invars <- names(fit1.glm$coef)
  check <- invars %w/o% row.names(out)
  if(length(check) != 0) cat(c(check, "removed due to singularities"), "\n")

  out.list <- list(formula1=formula1, formula2=formula2, results=out,empiricalResults=allResults, summary.coefs=summary.coefs,ANOD=anovfull.out,logLik=likelihood.out, WALD=WALD.out, aic=colMeans(aic.dat), aicEmpirical=aic.dat, effect=Effect)
  class(out.list) <- "hapBin"
  return(out.list)
}

