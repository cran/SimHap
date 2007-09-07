SNP2Haplo <- function(geno){
      
      numSNPs <- ncol(geno)-1
      newdat <- as.data.frame(matrix(nrow=nrow(geno), ncol=(2*numSNPs)+1))
      newdat[,1] <- geno[,1]
      newdat[,seq(from=2, to=ncol(newdat),by=2)] <- geno[,2:(numSNPs+1)]
      #print(newdat)
      
      for(j in seq(from=2, to=ncol(newdat),by=2)){

             newdat[,j] <- as.vector(newdat[,j], mode="character")        
             newdat[,j] <- replace(newdat[,j], newdat[,j]=="", NA)
             newdat[,j] <- replace(newdat[,j], substr(newdat[,j],1,1)==" ", NA)
             newdat[,j+1] <- substr(newdat[,j], nchar(newdat[,j])/2 + 1, nchar(newdat[,j]))
	     newdat[,j+1] <- replace(newdat[,j+1], newdat[,j+1]=="", NA)
             newdat[,j+1] <- replace(newdat[,j+1], substr(newdat[,j+1],1,1)==" ", NA)
	     
             newdat[,j] <- substr(newdat[,j], 1, nchar(newdat[,j])/2)
             newdat[,j] <- replace(newdat[,j], newdat[,j]=="", NA)
             newdat[,j] <- replace(newdat[,j], substr(newdat[,j],1,1)==" ", NA)

      }
      
      names(newdat)[1] <- names(geno)[1]
      names(newdat)[seq(from=2, to=ncol(newdat),by=2)] <- paste(names(geno)[2:(numSNPs+1)], "_1", sep="")
      names(newdat)[seq(from=2, to=ncol(newdat),by=2)+1] <- paste(names(geno)[2:(numSNPs+1)], "_2", sep="")
      
   return(newdat)
   }
   