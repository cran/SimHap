`organise.haplos` <-
function(geno) {

for(i in 1:nrow(geno)){
    for(j in 1:ncol(geno)){
      if(is.na(geno[i,j])) geno[i,j] <- c(" ")
    }
  }


  numSNPs <- ncol(geno)/2

# force columns of geno to have mode "character" to eliminate confusion
# with numeric alleles

  geno<-as.matrix(geno)
  for(i in 1:ncol(geno)){
     mode(geno[,i]) <- "character"
     }

  for (i in seq(length=numSNPs, from=1, by=2)){
       if (length(na.omit(unique(c(geno[,i],geno[,i+1])))) > 3)
           {
           stop("Length of allele cannot be greater than 2 characters.\n")
           }
       }

  # if SNP data are given as allele letters, create all.code,
  # a vector of the more frequent allele in each SNP
  # and encode geno into 0's and 1's

  geno2 <- geno
  for(i in 1:nrow(geno2)){
    for(j in 1:ncol(geno2)){
      if(geno2[i,j]== " " || geno2[i,j]== "") geno2[i,j] <- NA
    }
  }

      all.code <- vector("character", 0)  # most frequent allele in SNP
      all.code1 <- vector("character", 0) # least frequent allele in SNP
      for (i in seq(length=numSNPs, from=1, by=2)){
          # monomorphic SNPs
          if(length(unique(na.omit(c(geno2[,i],geno2[,i+1]))))==1){
            allele1 <- unique(na.omit(c(geno2[,i],geno2[,i+1])))[1]
            all.code <- c(all.code, allele1)
            all.code1 <- c(all.code1, allele1)

            geno[geno[,i]==allele1,i] <- 0
            geno[geno[,i+1]==allele1,i+1] <- 0
          }
          else {
            allele1 <- unique(na.omit(c(geno2[,i],geno2[,i+1])))[1]
            allele2 <- unique(na.omit(c(geno2[,i],geno2[,i+1])))[2]

            if (sum(na.omit(geno2[,i:(i+1)])==allele1) >= sum(na.omit(geno2[,i:(i+1)])==allele2)) {
              all.code <- c(all.code, allele1); all.code1<- c(all.code1, allele2)
            }
            else {
              all.code <- c(all.code, allele2); all.code1<- c(all.code1, allele1)
            }


            geno[geno[,i]==tail(all.code, n=1),i]<-0
            geno[geno[,i+1]==tail(all.code, n=1),i+1]<-0
            geno[geno[,i]==tail(all.code1, n=1),i]<-1
            geno[geno[,i+1]==tail(all.code1, n=1),i+1]<-1
          }

          #geno[!geno[,i]==0,i]<-1
          #geno[!geno[,i+1]==0,i+1]<-1
      }

#----------------------------------------------------------------------

  ID <- c(1:nrow(geno))
  haplo.labels <- haplo.label.numeric(0:(2^numSNPs-1),numSNPs=numSNPs)
  numHaplos <- nrow(haplo.labels)

  # Build data frames to hold haplotype data design

  wt <- NULL
  haploDM <- matrix(nrow=(nrow(geno)*2^(numSNPs-1)), ncol=numHaplos)
  DMindex <- 1
  hapMat <- matrix(nrow=(nrow(geno)*2^(numSNPs-1)), ncol=ncol(geno))
  hapMatIndex <- 1
  ID.vec <- matrix(nrow=(nrow(geno)*2^(numSNPs-1)),ncol=1)
  IDindex <- 1
  het.index <- rep(NA,numSNPs)


  # Construct design matrix
  # For phase ambiguous indivs add rows for each possible diplotype

  for(i in 1:nrow(geno)){
    for(j in 1:numSNPs){
       het.index[j] <- (.subset(geno,i,2*j-1)!=.subset(geno,i,2*j))
       }

    numHetero <- sum(het.index)

    # The rows of matrix myhaplos are possible diplotypes
    # for the ith subject.

    myhaplos <- make.haplo(geno[i,],het.index)
    numHaploComb <- nrow(myhaplos)

    #loop over haplo combos consistent w/ obs data
    for(j in 1:numHaploComb) {
      haploDM[DMindex,] <- codeDM(myhaplos[j,],haplo.labels)
      DMindex <- DMindex+1
    }


    for(j in 1:numHaploComb){
      ID.vec[IDindex] <- ID[i]
      IDindex <- IDindex+1
    }

    for(j in 1:nrow(myhaplos)) {
        hapMat[hapMatIndex,] <- myhaplos[j,]
        hapMatIndex <- hapMatIndex+1
    }

  }
  #end for loop over subjects

  haploDM <- haploDM[1:(DMindex-1),]
  hapMat <- hapMat[1:(hapMatIndex-1),]
  ID.vec <- ID.vec[1:(IDindex-1),]

  #Now renormalize the weights to sum to 1
  for(i in 1:(IDindex-1)) {
    wt[i] <- 1/sum(ID.vec==ID.vec[i])
  }

  # return the columns of haploDM with non-zero column sums

  myColSums <- colSums(haploDM)
  haploDM <- haploDM[,myColSums>0]

  haploDM <- data.frame(haploDM)

  n<-ncol(hapMat)
  hapMat2<-matrix(nrow=nrow(hapMat),ncol=2)

  hapMat2[,1] <- paste("h",hapMat[,1],sep="")
  hapMat2[,2] <- paste("h",hapMat[,n/2+1],sep="")
  for(i in 2:(n/2)) {
     hapMat2[,1] <- paste(hapMat2[,1],hapMat[,i],sep="")
     hapMat2[,2] <- paste(hapMat2[,2],hapMat[,i+n/2],sep="")
  }

  #Need to protect columns of hapMat2 from being coerced into factors
  #with the I() function.

  hapMat2 <- data.frame(haplo1=I(hapMat2[,1]),haplo2=I(hapMat2[,2]))

  hdmnames<-haplo.label(0:(2^numSNPs-1),numSNPs)

  # convert allele codes back to original labelling from input file
  # Allow conversion from single character numeric code to double-
  # character string.

  namevect <- vector()

    for (i in 1:length(hdmnames)) {
      namevect[i] <- ""
         for (j in 1:nchar(hdmnames[i])) {
            if(substr(hdmnames[i],j,j)==0) namevec <- all.code[j] else namevec <- all.code1[j]
            namevect[i] <- paste(namevect[i], namevec, sep="")
         }
    }
   hdmnames <- namevect

      # fix hapMat:
      HapNameVect <- vector()
      for (h in 1:length(hapMat2)) {
          for (i in 1:length(hapMat2[,h])) {
              HapNameVect[i] <- "h."

              for (j in 2:nchar(hapMat2[i,h])) {
                  if(substr(hapMat2[i,h],j,j)==0) HapNameVec <- all.code[j-1] else HapNameVec <- all.code1[j-1]

                  HapNameVect[i] <- paste(HapNameVect[i], HapNameVec, sep="")
              }
                hapMat2[i,h] <- HapNameVect[i]

          }

      }


  names(haploDM)<-paste("h", hdmnames[myColSums>0],sep=".")

  for(i in 1:nrow(hapMat2)){
    for(j in 1:ncol(hapMat2)){
       if(!(hapMat2[i,j]%in%names(haploDM))) hapMat2[i,j] <- "miss"
    }
  }
  return(list(haploDM=haploDM, hapMat=hapMat2, wt=wt, ID=ID.vec))

}

