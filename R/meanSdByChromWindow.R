
### names of sd.by.scan.chrom.window must be chromosome ids


meanSdByChromWindow <- function(sd.by.scan.chrom.window, sex)
{ 

  chromosomes <- names(sd.by.scan.chrom.window)
  n.chromosomes <- length(chromosomes)  
  
  bvals <- sd.by.scan.chrom.window[1:n.chromosomes]
  sds.chr <- vector("list", n.chromosomes)
  
  for(s in 1:n.chromosomes) 
  { 
  	if (chromosomes[s] != "X")
        {	  
		sds.chr[[s]] <- matrix(nrow=2, ncol=ncol(bvals[[s]]))
		rownames(sds.chr[[s]]) <- c("Mean", "SD")
		#
		for(i in 1:ncol(bvals[[s]])) 
		{ 
			sds.chr[[s]]["Mean",i] <- mean(bvals[[s]][,i], na.rm=TRUE)
			sds.chr[[s]]["SD",i] <- sd(bvals[[s]][,i], na.rm=TRUE)
		} 
        
        }
        else
        {
        	sds.chr[[s]] <- matrix(nrow=4, ncol=ncol(bvals[[s]]))
        	rownames(sds.chr[[s]]) <- c("Female Mean", "Male Mean", "Female SD", "Male SD")

        	# figure out which samples are M & which are F
        	fem.ind <- !is.na(sex) & sex == "F"  # these are the indices for the female samples
        	
        	
        	# loop through bins of X chr, take out females # MOVED fem.ind in the 4 rows within loop
        	for(i in 1:ncol(bvals[[s]])) 
        	{ 	
			 sds.chr[[s]]["Female Mean",i] <- mean(bvals[[s]][ fem.ind,i], na.rm=TRUE)
			 sds.chr[[s]]["Male Mean",i] <-   mean(bvals[[s]][!fem.ind,i], na.rm=TRUE)
			 sds.chr[[s]]["Female SD",i] <-   sd(  bvals[[s]][ fem.ind,i], na.rm=TRUE)
			 sds.chr[[s]]["Male SD",i] <-     sd(  bvals[[s]][!fem.ind,i], na.rm=TRUE)
		} 
        }  
  } 

  names(sds.chr)<- chromosomes

  return(sds.chr)
} 


