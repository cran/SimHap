`haplo.quant` <-
function(formula1, formula2, pheno, haplo, sim, effect="add", sub=NULL,predict_variable=NULL) {
  library(stats)
  call <- match.call()

  hapFreqs <- haplo$hapObject$final.freq
  haplo <- haplo$hapData
  
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
  
  #fit smaller, nested model (without haplotypes)
  fit2.lm <- eval(substitute(lm(formula2, data=pheno, subset=subset), list(subset=sub)))
  
  #log likelihood for smaller, nested model (without haplotypes)
  lnLsmall <- logLik(fit2.lm)
  
  for(i in 1:sim) {

    # determine weight for each individual and populate hapXs vectors
    for(j in 1:num_indivs) {

      # save processing time ... ;)
      current_numw <- num_weights[j,1]

      # invalid data case, weights for an indiv do not add up to one
      if(current_numw == 0) {
        print("NOTE: Intermediate result returned")
        print("Error. Weights do not sum to 1, indiv ID:")
        print(j)
        print(i)
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

  # perform loop through all iterations, constructing the dataframe and applying lm to each one
#------------------------------------
# initialise output

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
  rsquared.dat <- NULL
  likelihood.dat <- NULL
  predictedfit.dat <- NULL
  predictedse.dat <- NULL
  predictedMat.dat <- NULL
  aic.dat <- NULL
  lr.dat <- NULL
  lrt.dat <- NULL
  lnLbig.dat <- NULL
  
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
   if(effect=="rec") {
    # Recessive model
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
    
    # perform the lm with the current dataframe
    fit1.lm <- eval(substitute(lm(formula1, data=dataframe, subset=subset), list(subset=sub)))

    if(is.null(predict_variable)) {
    
         predicted <- "No variable was chosen to generate marginal means"
         predicted.mat <- NULL
         predicted.fit <- NULL
         predicted.se <- NULL
    }

    else {
  
         nls <- length(levs <- levels(as.factor(dataframe[, predict_variable])))
         tmp1 <- mean(dataframe[formula1_terms], na.rm=TRUE)
         tmp2 <- data.frame(matrix(rep(tmp1, nls), nrow=nls, byrow=TRUE))
         names(tmp2) <- names(tmp1)
         tmp2[, predict_variable] <- as.factor(levs)
	 dataframe[, predict_variable] <- as.factor(dataframe[, predict_variable])
	 formula.pred <- formula(paste(rownames(attr(terms(formula1),"factors"))[1],"~", paste(formula1_terms[1:length(formula1_terms)],collapse="+")))
  	 fit.pred <- eval(substitute(lm(formula=formula.pred, data=dataframe, subset=subset),list(subset=sub)))

         if(attr(fit1.lm$terms,"dataClasses")[grep(predict_variable, names(attr(fit1.lm$terms,"dataClasses")))]=="factor") {
               predicted.values <- predict(fit.pred, newdata=tmp2, se.fit=T)
               predicted.fit <- predicted.values$fit
               predicted.se <- predicted.values$se.fit
               predicted.mat <- NULL

               for(j in 1:length(predicted.values$fit)){
                        predicted.mat <- cbind(predicted.mat,predicted.values$fit[j], predicted.values$se.fit[j])
                        }
               } 
        else { 
               predicted <- "SNP variable was fitted with class 'numeric' but class 'factor' was required"
               predicted.mat <- NULL  
               predicted.fit <- NULL
               predicted.se <- NULL
               }
    } 
    
    # lm
    fit.lm <- as.data.frame(summary(fit1.lm)$coefficients)
    anov <- as.data.frame(anova(fit2.lm, fit1.lm))
    anovfull <- as.data.frame(anova(fit1.lm, test="Chisq"))

    #extract log-likelihood of model with haplotypes
    #likelihood <- logLik(fit1.lm)
    lnLbig <- logLik(fit1.lm)
    lnLbig.dat <- rbind(lnLbig.dat, lnLbig)
    lr <- -2*(lnLsmall[1]-lnLbig[1])
    lr.dat <- rbind(lr.dat, lr)
    lr.df <- attr(lnLbig, "df")-attr(lnLsmall,"df")
    lrt <- pchisq(lr,df=lr.df)
    lrt.dat <- rbind(lrt.dat, lrt)

    #likelihood.dat <- rbind(likelihood.dat, likelihood)
    aic <- AIC(fit1.lm)
    aic.dat <- rbind(aic.dat, aic)
    rsquared <- (summary(fit1.lm))$adj.r.squared
    rsquared.dat <- rbind(rsquared.dat, rsquared)
    predictedfit.dat <- rbind(predictedfit.dat,predicted.fit)
    predictedse.dat <- rbind(predictedse.dat,predicted.se)
    predictedMat.dat <- rbind(predictedMat.dat,predicted.mat)
    
    # add this row to anovfull.dat
    anovfullp.row <- anovfull$"Pr(>F)"
    anovfullp.dat <- rbind(anovfullp.dat, anovfullp.row)
    anovfulldf.row <- anovfull$Df
    anovfulldf.dat <- rbind(anovfulldf.dat, anovfulldf.row)

    # add this row to anov.dat
    anovp.row <- anov$"Pr(>F)"[2]
    anovp.dat <- rbind(anovp.dat, anovp.row)
    anovdf.row <- anov$Df
    anovdf.dat <- rbind(anovdf.dat, anovdf.row)
    anovresdf.row <- anov$"Res.Df"
    anovresdf.dat <- rbind(anovresdf.dat, anovresdf.row)

    # add this row to coef.dat
    coef.row <- fit.lm$Estimate
    #coef.dat <- rbind(coef.dat, coef.row)
    coef.dat <- rbind(coef.dat, coef.row)
    

    # extract some elements from the lm summary method and add row to p.dat
    pvals <- t(fit.lm[,ncol(fit.lm)])
    p.dat <- rbind(p.dat, pvals)

    # add a row to stderror.dat
    stderrors <- fit.lm$"Std. Error"
    stderror.dat <- rbind(stderror.dat, stderrors)

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

  row.names(coef.dat) <- NULL
  row.names(p.dat) <- NULL
  row.names(stderror.dat) <- NULL
  row.names(predictedfit.dat) <- NULL
  row.names(predictedse.dat) <- NULL
  row.names(predictedMat.dat) <- NULL
  row.names(aic.dat) <- NULL
  row.names(lr.dat) <- NULL
  row.names(lrt.dat) <- NULL
  row.names(lnLbig.dat) <- NULL

  coef.dat <- as.data.frame(coef.dat)
  p.dat <- as.data.frame(p.dat)
  stderror.dat <- as.data.frame(stderror.dat)
  anovdf.dat <- as.data.frame(anovdf.dat)
  anovp.dat <- as.data.frame(anovp.dat)

  anovfulldf.dat <- as.data.frame(anovfulldf.dat)
  anovfullp.dat <- as.data.frame(anovfullp.dat)
  predictedfit.dat <- as.data.frame(predictedfit.dat)
  predictedse.dat <- as.data.frame(predictedse.dat)
  
  if(is.null(predict_variable)) pred.out <- "No variable was chosen to generate marginal means"
  
  if(!is.null(predict_variable) && attr(fit1.lm$terms,"dataClasses")[grep(predict_variable, names(attr(fit1.lm$terms,"dataClasses")))]!="factor") pred.out <- "SNP variable was fitted with class 'numeric' but class 'factor' was required"
  
  if(!is.null(predictedMat.dat)) {
  
       predictedMat.dat <- as.data.frame(predictedMat.dat)
  
       if(effect=="add") {
           predicted <- round(mean(predictedMat.dat), digits=max(3, getOption("digits") - 3))
	   predvals <- rbind(predicted[1], predicted[3], predicted[5])
	   stderrs <- rbind(predicted[2], predicted[4], predicted[6])
	   pred.out <- as.data.frame(cbind(predvals, stderrs))
	   names(pred.out) <- c("Predicted Value", "Std. Error")
	   rownames(pred.out) <- c("0 copies", "1 copy", "2 copies")
           names(predictedMat.dat)[seq(2,ncol(predictedMat.dat),by=2)] <- "std.error"
           names(predictedMat.dat)[1] <- "0 copies"
   names(predictedMat.dat)[3] <- "1 copy"
   names(predictedMat.dat)[5] <- "2 copies"
           }
       if(effect=="dom") {
           predicted <- round(mean(predictedMat.dat), digits=max(3, getOption("digits") - 3))
           predvals <- rbind(predicted[1], predicted[3])
	   stderrs <- rbind(predicted[2], predicted[4])
	   pred.out <- as.data.frame(cbind(predvals, stderrs))
	   names(pred.out) <- c("Predicted Value", "Std. Error")
	   rownames(pred.out) <- c("0 copies", "1 or 2 copies")
           names(predictedMat.dat)[seq(2,ncol(predictedMat.dat),by=2)] <- "std.error"
   names(predictedMat.dat)[1] <- "0 copies"
   names(predictedMat.dat)[3] <- "1 or 2 copies"
   }
       if(effect=="rec") {
           predicted <- round(mean(predictedMat.dat), digits=max(3, getOption("digits") - 3))
           predvals <- rbind(predicted[1], predicted[3])
	   stderrs <- rbind(predicted[2], predicted[4])
	   pred.out <- as.data.frame(cbind(predvals, stderrs))
	   names(pred.out) <- c("Predicted Value", "Std. Error")
	   rownames(pred.out) <- c("0 or 1 copies", "2 copies")
           names(predictedMat.dat)[seq(2,ncol(predictedMat.dat),by=2)] <- "std.error"
   names(predictedMat.dat)[1] <- "0 or 1 copies"
   names(predictedMat.dat)[3] <- "2 copies"
   }
       
       } 
       
       
       
       
  #if(!is.null(predictedMat.dat)) names(predictedMat.dat)[seq(1,ncol(predictedMat.dat),by=2)] <- paste(seq(0,ncol(predictedMat.dat)/2-1,by=1),"copy")
  aic.dat <- as.data.frame(aic.dat)
  lr.dat <- as.data.frame(lr.dat)
  lrt.dat <- as.data.frame(lrt.dat)
  lnLbig.dat <- as.data.frame(lnLbig.dat)

  names(p.dat) <- names(coef.dat)
  names(stderror.dat) <- names(coef.dat)
  names(aic.dat) <- c("AIC")

  allResults <- list(Coef=coef.dat, Std.Error=stderror.dat, P.Value=p.dat)
  names(allResults$Coef) <- row.names(fit.lm)
  names(allResults$Std.Error) <- row.names(fit.lm)
  names(allResults$P.Value) <- row.names(fit.lm)

  sum.of.squares <- NULL
  for(i in 1:ncol(allResults$Std.Error)){
       sum.of.squares <- cbind(sum.of.squares,sum(allResults$Std.Error[,i]^2))
  }
  sum.of.squares <- as.data.frame(sum.of.squares)
  names(sum.of.squares) <- names(allResults$Std.Error)
  se1 <- sqrt(sum.of.squares/nrow(allResults$Std.Error))
  se2 <- sd(allResults$Coef)
  se.adj <- sqrt(se1^2 + se2^2)
  
  out.coef <- as.numeric(formatC(mean(coef.dat)))
  out.pval <- as.numeric(formatC(mean(p.dat)))
  out.se <- as.numeric(formatC(mean(stderror.dat)))
  #if(!is.null(predicted.dat)) predicted.vals <- formatC(mean(predicted.dat))
  
  summary.coefs <- data.frame(cbind(mean(allResults$Coef), as.numeric(as.vector(se.adj)), mean(allResults$P.Value)))
  names(summary.coefs) <- c("Coefficient", "Std.error", "P.Value")
  LRT.out1 <- cbind(round(mean(lnLbig.dat[,1]),digits=4), attr(lnLbig, "df"), round(mean(lr.dat),digits=4), signif((1-mean(lrt.dat)), digits=4))
  LRT.out2 <- cbind(round(lnLsmall[1],digits=4), attr(lnLsmall, "df"), "", "")
  row.names(LRT.out1) <- c("Full model")
  row.names(LRT.out2) <- c("Non-genetic")
  LRT <- rbind(LRT.out1, LRT.out2)
  LRT <- as.data.frame(LRT)
  names(LRT) <- c("logLik", "df", "LR", "P.Value")

  anovfull.out <- cbind(mean(anovfulldf.dat), formatC(mean(anovfullp.dat)))
  row.names(anovfull.out) <- row.names(anovfull)

  anovfull.out <- as.data.frame(anovfull.out)
  names(anovfull.out) <- c("DF", "P-Value")

  anov.out1 <- cbind(mean(anovresdf.dat[1]), "", "")
  row.names(anov.out1) <- c("1")
  anov.out2 <- cbind(mean(anovresdf.dat[2]), mean(anovdf.dat[2]), signif(mean(anovp.dat), digits=3))
  row.names(anov.out2) <- c("2")
  likelihood.out <- paste("'log Lik'", round(mean(lnLbig.dat), digits=3), paste("(df=", attr(lnLbig, "df"), ")", sep=""))
  rsquared.out <- rbind(round((summary(fit2.lm)$adj.r.squared),digits=3),round(mean(rsquared.dat),digits=3))
  rsquared.out <- as.data.frame(rsquared.out)
  names(rsquared.out) <- "Adjusted R-Squared"
  row.names(rsquared.out) <- c("Without Haplotypes", "Including Haplotypes")

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
  row.names(out) <- row.names(fit.lm)

  anov.out <- rbind(anov.out1, anov.out2)
  anov.out <- as.data.frame(anov.out)
  names(anov.out) <- c("Residual DF", "DF", "P-Value")

  if(effect=="add") Effect <- ("ADDITIVE")
  if(effect=="dom") Effect <- ("DOMINANT")
  if(effect=="rec") Effect <- ("RECESSIVE")


  "%w/o%" <- function(x,y) x[!x %in% y]
  invars <- names(fit1.lm$coef)
  check <- invars %w/o% row.names(out)
  if(length(check) != 0) cat(c(check, "removed due to singularities"), "\n")
  
  if(!is.null(predictedMat.dat)) {
     #predicted <- round(mean(predictedMat.dat), digits=max(3, getOption("digits") - 3))
     
     out.list <- list(formula1=formula1, formula2=formula2, results=out,empiricalResults=allResults, summary.coefs=summary.coefs, rsquared=rsquared.out,ANOD=anovfull.out,logLik=likelihood.out, LRT=LRT,predicted=pred.out, empiricalPredicted=predictedMat.dat, aic=mean(aic.dat), aicEmpirical=aic.dat, effect=Effect)
  } else {
         out.list <- list(formula1=formula1, formula2=formula2, results=out,empiricalResults=allResults, summary.coefs=summary.coefs, rsquared=rsquared.out,ANOD=anovfull.out,logLik=likelihood.out, LRT=LRT,predicted=pred.out, empiricalPredicted=predictedMat.dat, aic=mean(aic.dat), aicEmpirical=aic.dat, effect=Effect)
    }    
  class(out.list) <- "hapQuant"
  return(out.list)
}

