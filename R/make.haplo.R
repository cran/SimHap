`make.haplo` <-
function(geno.row,het.index){

    numSNPs <- length(het.index)
    num.hetero <- sum(na.omit(het.index))

    if(num.hetero<=1) {
        haplo <- matrix(nrow=1,ncol=length(geno.row))
        mid <- length(geno.row)/2
        haplo[1,(1:mid)] <- geno.row[2*(1:numSNPs)-1]
        haplo[1,(mid+1):length(geno.row)] <- geno.row[2*(1:numSNPs)]
    }
    else {
        is.hetero1 <- haplo.label.numeric(0:(2^(num.hetero-1)-1),numSNPs=num.hetero)
        is.hetero2 <- 1-is.hetero1
        haplo <- matrix(NA,ncol=2*numSNPs,nrow=2^(num.hetero-1))
        hap1 <- rep(FALSE,2*numSNPs)
        hap2 <- rep(FALSE,2*numSNPs)
        hap1[1:numSNPs] <- het.index
        hap2[(numSNPs+1):(2*numSNPs)] <- het.index

        for(i in 1:length(hap1)){
        if(is.na(hap1[i])) hap1[i] <- FALSE}

        haplo[,hap1] <- is.hetero1

        for(i in 1:length(hap2)){
        if(is.na(hap2[i])) hap2[i] <- FALSE}

        haplo[,hap2] <- is.hetero2

           if((numSNPs-num.hetero)>0) {
              is.homo <- geno.row[2*(1:numSNPs)][!het.index]
              is.homo <- matrix(rep(is.homo,2^(num.hetero-1)),ncol=(numSNPs-num.hetero),byrow=TRUE)
              hap1 <- rep(FALSE,2*numSNPs)             
              hap2 <- rep(FALSE,2*numSNPs)
              hap1[1:numSNPs] <- !het.index
              hap2[(numSNPs+1):(2*numSNPs)] <- !het.index

        for(i in 1:length(hap1)){
        if(is.na(hap1[i])) hap1[i] <- FALSE}

              haplo[,hap1] <- is.homo

        for(i in 1:length(hap2)){
        if(is.na(hap2[i])) hap2[i] <- FALSE}
              haplo[,hap2] <- is.homo

           }
        }

    return(haplo)
}

