
simulateGenotypeMatrix <- function(n.snps=10, 
                                   n.chromosomes=10, 
                                   n.samples=1000, 
                                   filename,
                                   file.type=c("gds", "ncdf"),
                                   silent = TRUE) {

    file.type <- match.arg(file.type)
    
    m <- n.snps*n.chromosomes # num of snps;  num of chromosomes*num of snps per chromosome
    n <- n.samples # this is the number of samples
    
    geno <- matrix(NA, m, n)
    
    ## simulate data to be put in the variables
    
    ## get allelic freq by random sampling of a uniform dist on (0,1)
    afreq <- runif(m)
    ## simulate genotypes by random sampling 2 gametes for each sample from binomial dist
    for (i in 1:m) geno[i,] <- rbinom(n, 2, afreq[i])
    
    ## simulate the missing calls for the genotypes
    mrate <- 0.05 # assumed to be the same across all snps
    nmiss <- rbinom(1, n*m, mrate) # find the number of missing over entire matrix
    imiss <- sample(1:(n*m),nmiss) # vectorized indices for calls to become missing
    geno[imiss] <- -1 # set the indices of particular genotypes to be missing
                                        # -1 is missing value
    table(geno)
    
    ## missing call rate for each snp
    miss_rate <- apply(geno, 1, function(x) length(x[x==-1])/length(x))
    summary(miss_rate)
    
    ## generate snp and sample integer ids (consecutive integer values)
    snpID <- 1:m 
    sampleID <- 1:n
    
    ## snp chromosome and position within chromosome values
    chromosome <- sort(rep(1:n.chromosomes, n.snps))
    position <- rep(seq(1:n.snps),n.chromosomes)
    
    ## create the file
    snp.annotation <- data.frame(snpID, chromosome, position)
    if (file.type == "gds") {
        var.data <- list(sample.id=sampleID, genotype=geno)
        nc <- .createGds(snp.annotation, filename, variables="genotype",
                         var.data=var.data)
    } else if (file.type == "ncdf") {
        var.data <- list(sampleID=sampleID, genotype=geno)
        nc <- .createNcdf(snp.annotation, filename, variables="genotype",
                          n.samples=length(sampleID), var.data=var.data)
    }
    .close(nc)
    
    if (!silent)
        return(table(geno)) 
    
}


