`infer.haplos` <-
function(geno){

  geno.ids <- geno[,1]
  geno <- geno[,-1]
  thresh=1/(2*nrow(geno)*10)
  numSNPs <- ncol(geno)/2
  haplo.object <- organise.haplos(geno)
  haplotest <- FALSE; ID.check <- rep(FALSE,length(haplo.object$ID))
  newhaploDM <- haplo.object$haploDM
  newhaploMat <- haplo.object$hapMat
  newID <- haplo.object$ID

  #Initial EM algorithm
  EM.results <- initEM(haplo.object)
  post.prob <- EM.results$wts
  initFreq <- EM.results$freq[!EM.results$freq<thresh] #haplos with non-zero frequency
  if(any(names(initFreq)=="miss")) initFreq <- initFreq[-which(names(initFreq)=="miss")]
  haplos.names <- names(newhaploDM)
  #haplos.names <- names(initFreq)
  zeroFreqHaplos <- names(EM.results$freq[EM.results$freq<thresh])
  if(sum(EM.results$freq<thresh)>0) { #then non-existent haplos need to be removed
    haplotest <- TRUE
    newhaploDM <- newhaploDM[,haplos.names]

    #We only want rows that sum to two - others must have involved
    #haplotypes with estimated frequency of zero
    
    #finalMatInd <- (rowSums(newhaploDM) == 2)

    #newhaploDM <- newhaploDM[finalMatInd,]
    #newhaploMat <- newhaploMat[finalMatInd,]

    #post.prob <- post.prob[finalMatInd]
    #newID <- newID[finalMatInd]
 
    # Re-calculate weights
    uniqueIDs <- unique(newID)
    IDsum <- rep(0,length(uniqueIDs))
    for(i in 1:length(uniqueIDs)) {
      IDsum[i] <- sum(post.prob[newID==uniqueIDs[i]])
    }
    for(i in 1:length(post.prob)) {
      post.prob[i] <- signif(post.prob[i]/IDsum[uniqueIDs==newID[i]], digits=5)
    }
  }

  # replace temp IDs (1 to n) with actual IDs from data file

  j <- as.vector(unique(geno.ids))
  newID2 <- vector(length=0)
    for(i in 1:length(newID)){
      for(k in 1:length(j)){
        if(newID[i]== k) newID2[i] <- j[k]
      } 
    }
  hap.freq <- EMloop(haplo.object, results=EM.results)
  hapMat <- cbind(newID2, newhaploMat, post.prob)
  names(hapMat) <- c("ID", "haplo1", "haplo2", "post.prob")
  return(list(hapMat=hapMat, hap.freq=hap.freq, initFreq=initFreq))
  #return(list(newID=newID, ID=newID2, hapMat=hapMat, hap.freq=hap.freq, initFreq=initFreq))
}

