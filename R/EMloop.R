`EMloop` <-
function(haplo.object, results)  {
 # Get the results of the last iteration

 results$freq <- results$freq[results$freq>0.001]
 if(any(names(results$freq)=="miss")) results$freq <- results$freq[-which(names(results$freq)=="miss")]
 freqnames<-rownames(results$freq)
 freq<-as.vector(results$freq)
 num.freq<-length(freq)
 weights<-results$wt
 ID <- haplo.object$ID
 #missing<-(weights<1)

   hapMat<-haplo.object$hapMat
   # matrix of T/F with (i,j)th element true if hapMat[i,k]==freqnames[j]
   # coerces to 1/0 s first before adding
   Xgen<-outer(hapMat[,1],freqnames,"==")+outer(hapMat[,2],freqnames,"==")
   #Xgen.mis <- as.matrix(Xgen[missing,])
 
   # Design matrices
   ones.mat<-matrix(1,nrow=(num.freq-1), ncol=(num.freq-1))
   num.haplo<-vector(length=num.freq)
   for (i in 1:num.freq){
      num.haplo[i]<-sum(weights*Xgen[,i])}

  # accommodate for Inf values in freq.deriv2 by taking a scalar multiple
  LAMBDA <- 10^100 # adjust as necessary
  freq.deriv <- 1/(LAMBDA*freq)
  freq.deriv2 <- (1/(LAMBDA*freq))^2 

  # If Inf values still exist, paste warning
  for(i in 1:length(freq.deriv2)){
    if(freq.deriv2[i] == Inf) print("Haplotype frequencies being generated might be too small due to many multiple hereozygotes. Consider reducing the number of SNPs being entered into the haplotype.") 
  }

   # Complete data expected information
   MIc.cov <- diag(num.haplo[1:(num.freq-1)]*freq.deriv2[1:(num.freq-1)])+freq.deriv2[num.freq]*num.haplo[num.freq]*ones.mat
   #return(MIc.cov)
   MIc <- MIc.cov
   varEM <- solve(MIc, LINPACK=TRUE, tol=1e-17)/LAMBDA^2
   numfreqs<-length(freqnames)
   varfreqs<-varEM
   cvec<-rep(-1,numfreqs-1)
   varlast<-t(cvec)%*%varfreqs%*%cvec
   freq.table<-data.frame(cbind(freqnames, freq, c(sqrt(diag(varfreqs)),sqrt(varlast))))

   names(freq.table) <- c("Haplotype", "Frequency", "Std.Error")

   freq.table$Frequency <- as.vector(freq.table$Frequency)
   freq.table$"Std.Error" <- as.vector(freq.table$"Std.Error")
   freq.table$Haplotype <- as.vector(freq.table$Haplotype)
   freq.table$Frequency <- as.numeric(freq.table$Frequency)
   freq.table$"Std.Error" <- as.numeric(freq.table$"Std.Error")
  
  
   ord <- order(-(freq.table$Frequency))
   freq.table <- freq.table[c(ord),]
   row.names(freq.table) <- sort(as.numeric(row.names(freq.table)))
  
   # Only include haplotypes with freq > 1e-05
   freq.table <- freq.table[freq.table$Frequency>1e-05,]
   return(freq.table)


}

