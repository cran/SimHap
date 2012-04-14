`freqTest` <- 
function(terms, freqs, n, effect){
      ntests <- length(which(freqs$Haplotype%in%terms))
      est.nums <- vector(mode="numeric", length=ntests)
      for(i in 1:ntests){
            if(effect=="add") est.nums[i] <- n*freqs$Frequency[which(freqs$Haplotype%in%terms)[i]]*(1-freqs$Frequency[which(freqs$Haplotype%in%terms)[i]])*2
	    else if(effect=="dom") est.nums[i] <- n*(1-(freqs$Frequency[which(freqs$Haplotype%in%terms)[i]])^2)
	    else est.nums[i] <- n*((freqs$Frequency[which(freqs$Haplotype%in%terms)[i]])^2)
      }
      est.nums.out <- as.data.frame(cbind(freqs$Haplotype[which(freqs$Haplotype%in%terms)], round(est.nums, digits=2)))
      names(est.nums.out) <- c("Haplotype", "Estimated Numbers")
      
      if(any(est.nums<=10)) warning(paste("Estimated number of individuals with haplotype:", est.nums.out$Haplotype[which(est.nums<=10)], "less than 10.","\n", "Statistical power to detect an association is very low, and collinearity may occur resulting in unreliable results.\n"), call.=FALSE)
      est.nums.out <- as.data.frame(cbind(freqs$Haplotype[which(freqs$Haplotype%in%terms)], round(est.nums, digits=2)))
  
      names(est.nums.out) <- c("Haplotype", "Estimated Numbers")
      est.nums.out
}
