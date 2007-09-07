SNP2Geno <- function(geno, baseline){

        numSNPs <- ncol(geno)-1
        if(length(baseline) != numSNPs) stop("Number of baseline genotypes does not match number of SNPs in geno")
        
        newdat <- as.data.frame(matrix(nrow=nrow(geno), ncol=(3*numSNPs)+1))
        newdat[,1] <- geno[,1]
        newdat[,seq(from=2, to=(3*numSNPs)+1,by=3)] <- geno[,2:(numSNPs+1)]
        newdat[,seq(from=2, to=(3*numSNPs)+1,by=3)+1] <- geno[,2:(numSNPs+1)]
        newdat[,seq(from=2, to=(3*numSNPs)+1,by=3)+2] <- geno[,2:(numSNPs+1)]
        

        k <- 1
        for(j in seq(from=2, to=(3*numSNPs)+1,by=3)){

            newdat[,j] <- as.vector(newdat[,j], mode="character")
            newdat[,j] <- replace(newdat[,j], newdat[,j]=="", NA)
            newdat[,j] <- replace(newdat[,j], substr(newdat[,j],1,1) == " ", NA)
            newdat[,j] <- replace(newdat[,j], newdat[,j]==baseline[k], 0)
            newdat[,j] <- replace(newdat[,j], nchar(newdat[,j])>1 & substr(newdat[,j], 1,nchar(newdat[,j])/2)!=substr(newdat[,j], nchar(newdat[,j])/2 + 1,nchar(newdat[,j])), 1)
            newdat[,j] <- replace(newdat[,j], nchar(newdat[,j])>1 & newdat[,j]!=baseline[k] & substr(newdat[,j], 1,nchar(newdat[,j])/2)==substr(newdat[,j], nchar(newdat[,j])/2 + 1,nchar(newdat[,j])), 2)
	    newdat[,j] <- as.numeric(newdat[,j])

            #dominant columns
            newdat[,j+1] <- as.vector(newdat[,j+1], mode="character")
            newdat[,j+1] <- replace(newdat[,j+1], newdat[,j+1]=="", NA)
            newdat[,j+1] <- replace(newdat[,j+1], substr(newdat[,j+1],1,1)==" ", NA)
            newdat[,j+1] <- replace(newdat[,j+1], newdat[,j+1]==baseline[k], 0)
            newdat[,j+1] <- replace(newdat[,j+1], nchar(newdat[,j+1])>1 & substr(newdat[,j+1], 1,nchar(newdat[,j+1])/2)!=substr(newdat[,j+1], nchar(newdat[,j+1])/2 + 1,nchar(newdat[,j+1])), 1)
            newdat[,j+1] <- replace(newdat[,j+1], nchar(newdat[,j+1])>1 & newdat[,j+1]!=baseline[k] & substr(newdat[,j+1], 1,nchar(newdat[,j+1])/2)==substr(newdat[,j+1], nchar(newdat[,j+1])/2 + 1,nchar(newdat[,j+1])), 1)
            newdat[,j+1] <- as.numeric(newdat[,j+1])
	    
            #recessive columns
            newdat[,j+2] <- as.vector(newdat[,j+2], mode="character")
            newdat[,j+2] <- replace(newdat[,j+2], newdat[,j+2]=="", NA)
            newdat[,j+2] <- replace(newdat[,j+2], substr(newdat[,j+2],1,1)==" ", NA)
            newdat[,j+2] <- replace(newdat[,j+2], newdat[,j+2]==baseline[k], 0)
            newdat[,j+2] <- replace(newdat[,j+2], nchar(newdat[,j+2])>1 & substr(newdat[,j+2], 1,nchar(newdat[,j+2])/2)!=substr(newdat[,j+2], nchar(newdat[,j+2])/2 + 1,nchar(newdat[,j+2])), 0)
            newdat[,j+2] <- replace(newdat[,j+2], nchar(newdat[,j+2])>1 & newdat[,j+2]!=baseline[k] & substr(newdat[,j+2], 1,nchar(newdat[,j+2])/2)==substr(newdat[,j+2], nchar(newdat[,j+2])/2 + 1,nchar(newdat[,j+2])), 1)
            newdat[,j+2] <- as.numeric(newdat[,j+2])
	    
            k <- k+1 
        }
        names(newdat)[1] <- names(geno)[1]
        names(newdat)[seq(from=2, to=(3*numSNPs)+1, by=3)] <- paste(names(geno)[2:(numSNPs+1)], "_add", sep="")
        names(newdat)[seq(from=2, to=(3*numSNPs)+1, by=3)+1] <- paste(names(geno)[2:(numSNPs+1)], "_dom", sep="")
        names(newdat)[seq(from=2, to=(3*numSNPs)+1, by=3)+2] <- paste(names(geno)[2:(numSNPs+1)], "_rec", sep="")
        return(newdat)
}


