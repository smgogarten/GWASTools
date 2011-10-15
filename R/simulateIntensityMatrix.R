#########
# simulateIntensityMatrix.r
#########

simulateIntensityMatrix<- 
function(n.snps=10, 
		 n.chromosomes=10, 
		 n.samples=1000, 
		 ncdf.filename, 
		 silent = TRUE) 
{
	
	m <- n.snps*n.chromosomes # total number of snps: num of chromosomes*num of snps per chromosome
	n <- n.samples # this is the number of samples
	
	
	## Generate data to be stored in NCDF file
	quantx <- matrix(NA, m, n)
	quanty <- matrix(NA, m, n)
	qual <- matrix(NA, m, n)
	
	# get intensities by random sampling of a normal dist on (0,1)
	# simulate genotypes by random sampling X & Y intensities for each sample from binomial dist
	het <- rbinom(m, 2, prob=0.3)
	for (i in 1:m) {
	if(het[i]==0) { meanX <- 0 } # homozygote BB
	if(het[i]==1) { meanX <- 1 } # heterozygote
	if(het[i]==2) { meanX <- 2 } # homozygote AA
	quantx[i,] <- abs(rnorm(n, mean=meanX) )
	quanty[i,] <- abs(rnorm(n, mean=abs(meanX-2)))
	qual[i,] <- rep(runif(1), n)
	}
	
	
	# simulate the missing data for the intensities
	mrate <- 0.05 # assumed to be constant over all snps
	nmiss <- rbinom(2, n*m, mrate) # find the number of missing over entire matrix
	imiss <- sample(1:(n*m), nmiss[1]) # vectorized indices for calls to become missing
	quantx[imiss] <- -1 # set the indices to be missing to have value NA
	quanty[imiss] <- -1
	qual[imiss] <- -1
	length(imiss) # total number of missing intensity values
	
	# generate snp and sample integer ids (consecutive integer values)
	snpid <- 1:m 
	sampleid <- 1:n
	
	# snp chromosome and position within chromosome values
	chr <- sort(rep(1:n.chromosomes, n.snps))
	pos <- rep(seq(1:n.snps),n.chromosomes)
	
	
	
	## Create & Load NCDF file with data
	# to create a ncdf file of samples and snps, first make dimensions
	sample <- dim.def.ncdf("sample", "id", sampleid)
	snp <- dim.def.ncdf("snp","count", snpid)
	
	# now create variables
	sampleID <- var.def.ncdf("sampleID", "id", sample, 0, prec="integer")
	position <- var.def.ncdf("position", "id", snp, -1, prec="integer")
	chromosome <- var.def.ncdf("chromosome", "id", snp, -1, prec="integer")
	quality <- var.def.ncdf("quality", "value", list(snp, sample), -9999, prec="float")
	X <- var.def.ncdf("X", "value", list(snp,sample), -9999, prec="float")
	Y <- var.def.ncdf("Y", "value", list(snp,sample), -9999, prec="float")
		 
	# then create the netCDF file itself
	ncquant <- create.ncdf(ncdf.filename, list(sampleID, position, chromosome, quality, X, Y))
	
	# now populate the ncdf variables with the above generated values
	put.var.ncdf(ncquant, "sampleID", sampleid)
	put.var.ncdf(ncquant, "position", pos)
	put.var.ncdf(ncquant, "chromosome", chr)
	put.var.ncdf(ncquant, "quality", qual)
	put.var.ncdf(ncquant, "X", quantx)
	put.var.ncdf(ncquant, "Y", quanty)
	
	close.ncdf(ncquant)

	if (!silent)
		return(list(het = het, nmiss = nmiss))
	

}
