#########
# simulateGenotypeMatrix.r
#########

simulateGenotypeMatrix<- 
function(n.snps=10, 
		 n.chromosomes=10, 
		 n.samples=1000, 
		 ncdf.filename,
		 silent = TRUE) 
{

	m <- n.snps*n.chromosomes # num of snps;  num of chromosomes*num of snps per chromosome
	n <- n.samples # this is the number of samples
	
	geno <- matrix(NA, m, n)
	
	# simulate data to be put in the variables
	
	# get allelic freq by random sampling of a uniform dist on (0,1)
	afreq <- runif(m)
	# simulate genotypes by random sampling 2 gametes for each sample from binomial dist
	for (i in 1:m) geno[i,] <- rbinom(n, 2, afreq[i])
	
	# simulate the missing calls for the genotypes
	mrate <- 0.05 # assumed to be the same across all snps
	nmiss <- rbinom(1, n*m, mrate) # find the number of missing over entire matrix
	imiss <- sample(1:(n*m),nmiss) # vectorized indices for calls to become missing
	geno[imiss] <- -1 # set the indices of particular genotypes to be missing
			# -1 is missing value
	table(geno)
	
	# missing call rate for each snp
	miss_rate <- apply(geno, 1, function(x) length(x[x==-1])/length(x))
	summary(miss_rate)
	
	# generate snp and sample integer ids (consecutive integer values)
	snpid <- 1:m 
	sampleid <- 1:n
	
	# snp chromosome and position within chromosome values
	chr <- sort(rep(1:n.chromosomes, n.snps))
	pos <- rep(seq(1:n.snps),n.chromosomes)
	
	# to create a ncdf file of samples and snps, first make dimensions
	sample <- dim.def.ncdf("sample", "id", sampleid)
	snp <- dim.def.ncdf("snp","count", snpid)
	
	# now create variables
	sampleID <- var.def.ncdf("sampleID", "id", sample, 0, prec="integer")
	genotype <- var.def.ncdf("genotype", "count", list(snp,sample), -1, prec="byte")
	position <- var.def.ncdf("position", "count", snp, -1, prec="integer")
	chromosome <- var.def.ncdf("chromosome", "id", snp, -1, prec="integer")
	
	# then create the netCDF file itself
	ncgeno <- create.ncdf(ncdf.filename, list(sampleID,genotype, position, chromosome))
	
	# now populate the ncdf variables with the above generated values
	put.var.ncdf(ncgeno, "sampleID", sampleid)
	put.var.ncdf(ncgeno, "genotype", geno)
	put.var.ncdf(ncgeno, "position", pos)
	put.var.ncdf(ncgeno, "chromosome", chr)
	
	close.ncdf(ncgeno)
	
	if (!silent)
		return(table(geno)) 
	
}


