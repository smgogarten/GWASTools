
simulateIntensityMatrix <- function(n.snps=10, 
                                    n.chromosomes=10, 
                                    n.samples=1000, 
                                    filename, 
                                    file.type=c("gds", "ncdf"),
                                    silent = TRUE) {
    
    file.type <- match.arg(file.type)
    
    m <- n.snps*n.chromosomes # total number of snps: num of chromosomes*num of snps per chromosome
    n <- n.samples # this is the number of samples
    
    
    ## Generate data to be stored in NCDF file
    quantx <- matrix(NA, m, n)
    quanty <- matrix(NA, m, n)
    qual <- matrix(NA, m, n)
    
    ## get intensities by random sampling of a normal dist on (0,1)
    ## simulate genotypes by random sampling X & Y intensities for each sample from binomial dist
    het <- rbinom(m, 2, prob=0.3)
    for (i in 1:m) {
	if(het[i]==0) { meanX <- 0 } # homozygote BB
	if(het[i]==1) { meanX <- 1 } # heterozygote
	if(het[i]==2) { meanX <- 2 } # homozygote AA
	quantx[i,] <- abs(rnorm(n, mean=meanX) )
	quanty[i,] <- abs(rnorm(n, mean=abs(meanX-2)))
	qual[i,] <- rep(runif(1), n)
    }
    
    
    ## simulate the missing data for the intensities
    mrate <- 0.05 # assumed to be constant over all snps
    nmiss <- rbinom(2, n*m, mrate) # find the number of missing over entire matrix
    imiss <- sample(1:(n*m), nmiss[1]) # vectorized indices for calls to become missing
    quantx[imiss] <- -1 # set the indices to be missing to have value NA
    quanty[imiss] <- -1
    qual[imiss] <- -1
    length(imiss) # total number of missing intensity values
    
    ## generate snp and sample integer ids (consecutive integer values)
    snpID <- 1:m 
    sampleID <- 1:n
    
    ## snp chromosome and position within chromosome values
    chromosome <- sort(rep(1:n.chromosomes, n.snps))
    position <- rep(seq(1:n.snps),n.chromosomes)
    
    ## create the file
    snp.annotation <- data.frame(snpID, chromosome, position)
    if (file.type == "gds") {
        var.data <- list(sample.id=sampleID, quality=qual, X=quantx, Y=quanty)
        nc <- .createGds(snp.annotation, filename, 
                         variables=c("quality", "X", "Y"),
                         var.data=var.data)
    } else if (file.type == "ncdf") {
        var.data <- list(sampleID=sampleID, quality=qual, X=quantx, Y=quanty)
        nc <- .createNcdf(snp.annotation, filename,
                          variables=c("quality", "X", "Y"),
                          n.samples=length(sampleID), var.data=var.data)
    }   
    .close(nc)
    
    if (!silent)
        return(list(het = het, nmiss = nmiss))
    

}
