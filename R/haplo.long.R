`haplo.long` <- 
function(fixed, random, pheno, haplo, cor="corCAR1", value=0.2, form=~1, sim, effect="add", sub=NULL) {
  
  library(stats)
  call <- match.call()
  library(nlme)
  hapFreqs <- haplo$hapObject$final.freq
  haplo <- haplo$hapData
  
  formula_fixed <- formula(fixed)
  formula_random <- formula(random)
  
  formula_fixednofactors <- formula_fixed
  formula_fixedterms <- attr(terms(formula_fixednofactors), "term.labels")
  
  if(any(regexpr(":", formula_fixedterms)!=-1)){
        formula_fixedterms <- formula_fixedterms[-which(regexpr(":", formula_fixedterms)!=-1)]
  }
  
  
  if(any(regexpr("factor", formula_fixedterms)==1)) {
        formula_fixedterms[which(regexpr("factor", formula_fixedterms)==1)] <- substr(formula_fixedterms[which(regexpr("factor", formula_fixedterms)==1)],8,nchar(formula_fixedterms[which(regexpr("factor", formula_fixedterms)==1)])-1)
  }   
  #else formula_fixedterms <- attr(terms(formula_fixednofactors), "term.labels")
  
  freq.estnums <- freqTest(terms=formula_fixedterms, freqs=hapFreqs, n=length(unique(haplo[,1])), effect=effect)
  
  num_indivs <- length(unique(as.numeric(haplo[,1])))
  #num_indivs <- as.numeric(haplo[nrow(haplo),1])

  # first column retains number of non-zero probabilities for individual
  # second column holds the current iteration index for the next weight change
  num_weights <- matrix(0, nrow=num_indivs, ncol=2)

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
  print("* Constructing dataframes and performing linear mixed effects model for each simulation ...")
  # This section constructs the dataframe for each iteration in preparation for lme

  haplo_table <- table(c(haplo[,2],haplo[,3]))
  num_haplos <- dim(haplo_table)
  names_haplos <- names(haplo_table)
  

  # prepare the reusable dataframe container ... 
  dataframe_extra <- matrix(0, nrow=num_indivs, ncol=num_haplos)

  # perform loop through all iterations, constructing the dataframe and applying glm to each one
#------------------------------------
  coef.dat <- NULL
  p.dat <- NULL
  stderror.dat <- NULL
  out <- NULL
  output <- NULL
  pvals <- NULL
  stderrors <- NULL
  
  anov.out <- NULL
  anov.out1 <- NULL
  anov.out2 <- NULL
  anovfull.out <- NULL
  anovfullp.dat <- NULL
  anovfulldf.dat <- NULL
  
  aic.dat <- NULL
  #lr.dat <- NULL
  #lrt.dat <- NULL
  lnLbig.dat <- NULL
  # the dynamic point at which the algorithm reports progress
  five_percent <- round(sim*0.05)
  report <- five_percent

  # main loop

#############################################

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
      colnames(dataframe_extra_expand)[j] <- paste(names_haplos[j])}

    dataframe <- as.data.frame(cbind(pheno, dataframe_extra_expand))
    dataframe <- dataframe[complete.cases(dataframe[formula_fixedterms]),]


    # perform the lme model with the current dataframe
    # lme

    if(cor=="corCAR1") {
    if(length(sub)==0) {
      fit1.lme <- lme(fixed=fixed, data=dataframe, random=random, correlation=corCAR1(form= form), na.action=na.omit)
      }
    
    else {
      fit1.lme <- eval(substitute(lme(fixed=fixed, data=dataframe, random=random, correlation=corCAR1(value=value, form= form), 
      na.action=na.omit, subset=subset),list(subset=sub)))
      }
    }


    if(cor=="corAR1") {
    if(length(sub)==0) {
      fit1.lme <- lme(fixed=fixed, data=dataframe, random=random, correlation=corAR1(value=value, form= form), na.action=na.omit)
      }
    
    else {
      fit1.lme <- eval(substitute(lme(fixed=fixed, data=dataframe, random=random, correlation=corAR1(value=value, form= form), 
      na.action=na.omit, subset=subset),list(subset=sub)))
      }
    }

    if(cor=="corCompSymm") {
    if(length(sub)==0) {
      fit1.lme <- lme(fixed=fixed, data=dataframe, random=random, correlation=corCompSymm(value=value, form=form), na.action=na.omit)
      }
    
    else {
      fit1.lme <- eval(substitute(lme(fixed=fixed, data=dataframe, random=random, correlation=corCompSymm(value=value, form=form), 
      na.action=na.omit, subset=subset),list(subset=sub)))
      }
    }

    sum.lme <- as.data.frame(summary(fit1.lme)$tTable)
    anovfull <- as.data.frame(anova(fit1.lme, test="Chisq"))
    anovfullp.row <- anovfull$"p-value"
    anovfullp.dat <- rbind(anovfullp.dat, anovfullp.row)
    anovfulldf.row <- anovfull$numDF
    anovfulldf.dat <- rbind(anovfulldf.dat, anovfulldf.row)
  
    lnLbig <- logLik(fit1.lme)
    lnLbig.dat <- rbind(lnLbig.dat, lnLbig)

    aic <- AIC(fit1.lme)
    aic.dat <- rbind(aic.dat, aic)

    # add this row to coef.dat

      coef.row <- sum.lme$Value
      coef.dat <- rbind(coef.dat, coef.row)
      pvals <- sum.lme$"p-value"
      p.dat <- rbind(p.dat, pvals)


      stderrors <- sum.lme$Std.Error
      stderror.dat <- rbind(stderror.dat, stderrors)

    # extract some elements from the lme summary method and add row to p.dat
    # report on progress
    if(i==report) {
      percentage <- report * 100 / sim
      print(paste(percentage, "%"))
      report <- report + five_percent
    }
  } 

  # nullify the row names
  row.names(coef.dat) <- NULL
  row.names(p.dat) <- NULL
  row.names(stderror.dat) <- NULL
  row.names(aic.dat) <- NULL
  row.names(lnLbig.dat) <- NULL
  row.names(anovfulldf.dat) <- NULL
  row.names(anovfullp.dat) <- NULL
 
  aic.dat <- as.data.frame(aic.dat)
  lnLbig.dat <- as.data.frame(lnLbig.dat)

  coef.dat <- as.data.frame(coef.dat)
  p.dat <- as.data.frame(p.dat)
  stderror.dat <- as.data.frame(stderror.dat)
  anovfulldf.dat <- as.data.frame(anovfulldf.dat)
  anovfullp.dat <- as.data.frame(anovfullp.dat)
  
  names(coef.dat) <- row.names(sum.lme)
  names(p.dat) <- row.names(sum.lme)
  names(stderror.dat) <- row.names(sum.lme)
  names(aic.dat) <- c("AIC")
  
  allResults <- list(Coef=coef.dat, Std.Error=stderror.dat, P.Value=p.dat)
  names(allResults$Coef) <- row.names(sum.lme)
  names(allResults$Std.Error) <- row.names(sum.lme)
  names(allResults$P.Value) <- row.names(sum.lme)
 
  sum.of.squares <- NULL
  for(i in 1:ncol(allResults$Std.Error)){
       sum.of.squares <- cbind(sum.of.squares,sum(allResults$Std.Error[,i]^2))
  }
  sum.of.squares <- as.data.frame(sum.of.squares)
  names(sum.of.squares) <- names(allResults$Std.Error)
  se1 <- sqrt(sum.of.squares/nrow(allResults$Std.Error))
  se2 <- sd(allResults$Coef)
  se.adj <- sqrt(se1^2 + se2^2)
  
  output$coef <- formatC(mean(coef.dat))
  output$pval <- formatC(mean(p.dat))
  output$se <- formatC(mean(stderror.dat))
 
  summary.coefs <- data.frame(cbind(mean(allResults$Coef), as.numeric(as.vector(se.adj)), mean(allResults$P.Value)))
  names(summary.coefs) <- c("Coefficient", "Std.error", "P.Value")
  anovfull.out <- cbind(mean(anovfulldf.dat), formatC(mean(anovfullp.dat)))
  row.names(anovfull.out) <- row.names(anovfull)

  anovfull.out <- as.data.frame(anovfull.out)
  names(anovfull.out) <- c("DF", "P-Value")
  likelihood.out <- paste("'log Lik'", round(mean(lnLbig.dat), digits=3), paste("(df=", attr(lnLbig, "df"), ")", sep=""))
  
  print("  Done")
  # ****************************

  # Arrange the output data
  
  for(i in 1:ncol(coef.dat)){

    output$coef.CI[i] <- paste("(",formatC(quantile(coef.dat[,i], probs=c(0.025), na.rm=T)),",",formatC(quantile(coef.dat[,i], probs=c(0.975), na.rm=T)),")", sep="")
    
    output$pval.CI[i] <- paste("(",formatC(quantile(p.dat[,i], probs=c(0.025), na.rm=T)),",",formatC(quantile(p.dat[,i], probs=c(0.975), na.rm=T)),")", sep="")
    output$se.CI[i] <- paste("(",formatC(quantile(stderror.dat[,i], probs=c(0.025), na.rm=T)),",",formatC(quantile(stderror.dat[,i], probs=c(0.975), na.rm=T)),")", sep="")
    
  }

  out <- data.frame(cbind(output$coef, output$coef.CI, output$se, output$se.CI, output$pval, output$pval.CI))

  names(out) <- c("Coef", "Coef.quantiles", "Std.Error", "Std.Error.quantiles", "P.Val", "P.Val.quantiles")
  if(effect=="add") Effect <- ("ADDITIVE")
  if(effect=="dom") Effect <- ("DOMINANT")
  if(effect=="rec") Effect <- ("RECESSIVE")
  
  out.list <- list(fixed_formula=formula_fixed, random_formula=formula_random, results=out,empiricalResults=allResults, summary.coefs=summary.coefs, ANOD=anovfull.out,logLik=likelihood.out, AIC=mean(aic.dat), aicEmpirical=aic.dat, corStruct=cor, effect=Effect)
  class(out.list) <- "hapLong"
  return(out.list)

}


