
# function for sliding window approach to finding b allele freq outlier detection
# takes all snps for each sample and slide the window along the data as such

# v2: writes the sample nums to output sd.chr list of matrices
# baf.mean.sd.windows revised by Cathy Laurie from slide.baf.v2 by Caitlin McHugh: 
# adds use of all probes that are not intensity only, not failed (missing=100%) and not called as homozygote (leaving hets and sporadic missing)
# and adds mean of each window in addition to sd


# inputs:
#  nbins: a vector of length 23 w entries for the number of bins for each chromosome.  recommended is approximately 8 segments/chr; must be an even number
#  intenData: IntensityData object with b allele freq values
#  genoData: GenotypeData object
#   snp.exclude:  a vector of "int.id"s to be exclude (i.e. not intensity only and not 100% missing)
#   gcall:  which genotype calls to include in variance calculation

# output: 
#   two lists of matrices, one for SD and one for the mean of BAF
#   each element of each list corresponding to a chromosome, each matrix dims number samples x number windows
#   entry [[i]][j,k] is the standard dev for sample j, window k in the ith chromosome


sdByScanChromWindow <- function(
	 intenData, 
	 genoData=NULL,
         var="BAlleleFreq", nbins=NULL,
         snp.exclude=NULL, return.mean=FALSE,
         incl.miss=TRUE, incl.het=TRUE, incl.hom=FALSE) 
{	
        # check gcall
        if (is.null(genoData) & !(incl.miss & incl.het & incl.hom)) {
          stop("genoData must be provided if any of incl.miss, incl.het, incl.hom is FALSE")
        }

        # check that intenData has BAF
        if (!hasVariable(intenData, varname=var)) stop(paste(var, "not found in intenData"))
  
        # check that dimensions of intenData and genoData are equal
        intenSnpID <- getSnpID(intenData)
        intenScanID <- getScanID(intenData)
        if (!is.null(genoData)) {
          genoSnpID <- getSnpID(genoData)
          if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
          genoScanID <- getScanID(genoData)
          if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")
        }
  
	nsamp <- length(intenScanID)  #nsamp is number of samples
	#
	chroms <- getChromosome(intenData, char=TRUE)
	pos <- getPosition(intenData)
	#
	## select chromosomes to run (only from 1:23)
	chromosomes <- unique(chroms)
	chromosomes <- chromosomes[is.element(chromosomes, c(autosomeCode(intenData), "X"))]
	#
	n.chrom <- length(chromosomes)
	
        # check nbins
        if (is.null(nbins)) {
          nbins <- rep(2, n.chrom)
        } else {
          if (length(nbins) != n.chrom)
            stop("The length of nbins must be equal to one more than the number of autosomes in the genotype netcdf file")
        }
        chk <- rep(NA,n.chrom)
        for(i in 1:n.chrom) chk[i] <- nbins[i]%%2 == 0
        if(!all(chk)) stop("each element of nbins must be an even number")
	
        # logical vector to select snps that do not include "snp.exclude"
        inc <- !is.element(intenSnpID, snp.exclude)
        # index for all snps
        m <- length(intenSnpID)
        index <- 1:m

	# fill a list w element for each chromosome w matrices w two cols: start index, end index for snps in that chr
	# for each window
	inds.chr <- vector("list", n.chrom)
	
	# create list to store SD values for each chromosome
	sd.chr <- vector("list", n.chrom)
        if (return.mean) mn.chr <- vector("list", n.chrom)
	
        totalsnps <- 0
        start.ind <- 1

	# get start and count row indices for each chromosome type

	indices <- matrix(NA, n.chrom, 2, dimnames=list(chromosomes, c("start","stop")))
	for(i in 1:n.chrom) 
		indices[i,] <- range(which(is.element(chroms, chromosomes[i])))
	
	# get snps per chromosome
	spc <- apply(indices, 1, function(x) x[2]-x[1]+1)
	
	## Get the starting and ending positions for each bin
	
	for(s in 1:n.chrom) 
	{
		#if(nbins[s]%%2 != 0) { nbins[s] <- nbins[s]+1 } # ensure number of bins is an even number
	
		#current.chrom <- chromosomes[s]
		
		nsnps.chr <- spc[s] 		# this is the number of snps for chromosome s
		nsnp.bin <- nsnps.chr/nbins[s]  # this is the number of snps per bin (can be fractional)
		nwin <- nbins[s]-1 		# this is number of windows for chromosome s
	
		start.ind <- indices[s,1]
		
		inds.chr[[s]] <- matrix(nrow=nwin, ncol=2) # col 1holds start index, col 2 holds end index for snps
		for(w in 1:nwin) 
		{
			inds.chr[[s]][w,1] <- start.ind
			inds.chr[[s]][w,2] <- start.ind+2*nsnp.bin
			start.ind <- start.ind + nsnp.bin
		}
		if (inds.chr[[s]][nwin,2] != indices[s,2]) 
			inds.chr[[s]][nwin,2] <- indices[s,2] 
		
		sd.chr[[s]] <- matrix(NA, nrow=nsamp, ncol=nwin)
                if (return.mean) mn.chr[[s]] <- matrix(NA, nrow=nsamp, ncol=nwin)
	} 

	
	# now we have a list of matrices for each chr that lists the start & end index for snps for each window
	# and we have a list of matrices for each chr of dim number samples x number of windows to store the sd values
	

	# loop through all samples i in 1:n
	for(i in 1:nsamp) 
	{ 
		# get b allele freq for all snps for sample i
                bv <- getVariable(intenData, varname=var, snp=c(1,-1), scan=c(i,1))
		#
		# get genotype data for all snps for sample i
                if (!is.null(genoData)) {
                  gv <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))
                }
		#
		for (j in 1:n.chrom) 
		{	
			for (k in 1:nrow(inds.chr[[j]])) 
			{ 
				start.ind <- inds.chr[[j]][k,1]
				end.ind <-   inds.chr[[j]][k,2]
				#
				balleles <- bv[start.ind:end.ind]
                                inc.bin <- inc[start.ind:end.ind] # logical for what snps to include
				#
				# calculate SD for genotype calls in gcall only; store in matrix entry [[i]][j,k]
                                if (!is.null(genoData)) {
                                  genos <- gv[start.ind:end.ind]
                                  sel <- rep(FALSE, length(balleles))
                                  if (incl.miss) {
                                    sel <- sel | is.na(genos)
                                  }                      
                                  if (incl.het) {
                                    sel <- sel | (!is.na(genos) & genos == 1)
                                  }                       
                                  if (incl.hom) {
                                    sel <- sel | (!is.na(genos) & (genos == 0 | genos == 2))
                                  }
                                } else {
                                  sel <- rep(TRUE, length(balleles))
                                }
                                
				sd.chr[[j]][i,k]<- sd(balleles[sel & inc.bin], na.rm=TRUE)
                                if (return.mean) mn.chr[[j]][i,k]<- mean(balleles[sel & inc.bin], na.rm=TRUE) 
				
			} # end loop through windows
			
		}  # end loop through chromosomes

	}  # end loop through samples
        
	# add scan IDs to the values being returned
	for(s in 1:n.chrom) {
                rownames(sd.chr[[s]]) <- intenScanID
                if (return.mean) rownames(mn.chr[[s]]) <- intenScanID
        }
        names(sd.chr) <- chromosomes
        if (return.mean) names(mn.chr) <- chromosomes

        if (return.mean) {
          return(list("sd"=sd.chr, "mean"=mn.chr))
        } else {
          return(sd.chr)
        }
}




