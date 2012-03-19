#####
# Cox proportional hazards
#####


assocTestCPH <- function(
	genoData,	# GenotypeData object containing sex and phenotypes
	event,	# name of variable in genoData for event to analyze
	time.to.event,		# name of variable in genoData for time to event
	covars,		# vector of covariate terms for model (can include interactions as 'a:b', main effects correspond to variable names in genoData)
	factor.covars,		# vector of names of covariates to be converted to factor
	scan.chromosome.filter = NULL,  # matrix of T/F for scan by chromosome for chromosome-specific selection of samples
        scan.exclude = NULL,
	maf.filter = FALSE,  # whether to filter results returned using maf  > 75/2n where n = number of events
	GxE = NULL,     # name of the covariate to use for E if genotype-by-environment (i.e. SNP:E) model is to be analyzed, in addition to the main effects (E can be a covariate interaction)
	stratum = NULL,  # name of variable to stratify on for a stratified analysis (use NULL if no stratified analysis needed)
	chromosome.set = NULL, 	# vector of chromosome numbers (corresponding to format of "chromosome" in genoData - i.e. integer codes)	
	block.size = 5000,	# number of SNPs from a given chromosome to read in one block from genoData
        verbose = TRUE,
        outfile = NULL
){

# v2 captures warnings and does not save the result if there is one
# v3 adds an indicator argument for whether to use the MAF criterion before fitting the model
# v3 also adds capture of error in addition to warning
# v4 saves number of events for each SNP
# v5 adds option for including one genotype-by-environment interaction term & make missing  value return for genotypes general
# v6 change returned data from matrix to data.frame to avoid floating point rounding

	############# preliminaries  ############################

	cvnames <-  unique(unlist(strsplit(covars,"[*:]"))) # names of covariate variables
	if(!all(is.element(c(event,time.to.event,cvnames), getScanVariableNames(genoData)))) stop("required phenotypes are not in genoData")

	chrom <- getChromosome(genoData)
	snpID <- getSnpID(genoData)
        scanID <- getScanID(genoData)

        # chromosome.set check
        if(is.null(chromosome.set)){
          chromosome.set <- unique(chrom)    
          allchrom <- TRUE   
	} else {
          stopifnot(all(chromosome.set %in% unique(chrom)))
          allchrom <- FALSE
        }

	sample.dat <- getScanVariable(genoData, unique(c("sex", cvnames, event, time.to.event, stratum)))  # get necessary variables

	if(!all(is.element(sample.dat$sex,c("M","F")))) stop("sex must be coded as M and F")

        # set or check scan.chromosome.filter
        if (is.null(scan.chromosome.filter)) {
          scan.chromosome.filter <- matrix(TRUE, nrow=length(scanID), ncol=max(chrom),
                                           dimnames=list(scanID, 1:max(chrom)))
          scan.chromosome.filter[sample.dat$sex=="F", YchromCode(genoData)] <-  FALSE
        } else {
          if(length(scanID) != nrow(scan.chromosome.filter)) stop("number of rows of scan.chromosome.filter does not match number of scans in genoData")
          if(!all(rownames(scan.chromosome.filter)==scanID)) stop("scan.chromosome.filter rows are out of order")
          if(!all(is.element(colnames(scan.chromosome.filter), unique(chrom)))) stop("column names of scan.chromosome.filter must be chromosome numbers")
          if(any(scan.chromosome.filter[sample.dat$sex=="F", YchromCode(genoData)])) stop("scan.chromosome.filter must be all FALSE for Y chromosome for females")
        }
	# if sex is in the model, make scan.chromosome.filter=FALSE for all subjects for the Y chr
	if(is.element("sex", cvnames))  scan.chromosome.filter[,YchromCode(genoData)] <- FALSE

        if (!is.null(scan.exclude)) {
          scan.chromosome.filter[(scanID %in% scan.exclude),] <- FALSE
        }
       
	if(!is.null(GxE)) {
          if( length(GxE)!=1 | !is.element(GxE, covars)) stop("GxE should be a single covariate")
        }
	
	if(!is.null(stratum)) {
          if(length(stratum)>1) stop("length of stratum must equal one")
        }

	model <- as.formula(paste("surv ~ ", paste(covars, collapse=" + "), " + gtype"))
	if(!is.null(GxE)) model2 <- as.formula(paste("surv ~ ", paste(covars, collapse=" + "), " + gtype +", GxE, ":gtype", sep=""))
	if(!is.null(stratum)) {
		model <- as.formula(paste("surv ~ ", paste(covars, collapse=" + "), " + gtype + strata(", stratum, ")", sep=""))
		if(!is.null(GxE)) model2 <- as.formula(paste("surv ~ ", paste(covars, collapse=" + "), " + gtype +", GxE, ":gtype+ strata(",stratum,")", sep=""))
	}

	# make matrix of chrom start and count
        c <- unique(chrom)
        nchr <- length(c)
        dup <- duplicated(chrom)
        x <- table(chrom)
        chrom.info <- matrix(NA, nchr,3)
        dimnames(chrom.info)[[2]] <- c("chrom", "start", "count")
        chrom.info[,"chrom"] <- c
        chrom.info[,"start"] <- which(!dup)
        chrom.info[,"count"] <- as.vector(x)

	# make data.frames for storing results
	# For main effects model (no SNP interaction)
        nsnp <- length(chrom[is.element(chrom, chromosome.set)])
        tmp <- rep(NA,nsnp)
        res <- data.frame(index=tmp, snpID=tmp, chr=tmp, maf=tmp, mafx=tmp, beta=tmp, se=tmp, z=tmp, pval=tmp, warned=tmp, n.events=tmp)
	# for model with SNP interaction
	if(!is.null(GxE)){
		res2 <- data.frame(index=tmp, snpID=tmp, chr=tmp, maf=tmp, mafx=tmp, warned=tmp, n.events=tmp, ge.lrtest=tmp, ge.pval=tmp)
		} else { 
		res2 <- NULL
	 }

	# make factors where necessary
        y <- length(factor.covars)
        for(q in 1:y) sample.dat[,factor.covars[q]] <- factor(sample.dat[,factor.covars[q]])


        k <- 1  # counter for row in res

	####################  End preliminaries  ##############################

	####################  Loop over chromosomes  #######################################

	for(chr in chromosome.set){
		
		sel <- scan.chromosome.filter[, chr]  # logical for selecting samples
		pheno <- sample.dat[sel,]  # select samples

		# get matrix of block info for a given chromosome (blk)

                cs <- chrom.info[chrom.info[,"chrom"] == chr, "start"]
                cn <- chrom.info[chrom.info[,"chrom"] == chr, "count"]
                m <- floor(cn/block.size)
                r <- cn - m*block.size

                if(m==0) nb <- 1
                if(m>0 & r<=1) nb <- m
                if(m>0 & r>1)  nb <- m + 1

                blk <- matrix(NA, nb, 3)
                dimnames(blk)[[2]] <- c("index", "start", "count")
                blk[,"index"] <- 1:nb
                blk[1,"start"] <- cs

                if(nb==1) { blk[1,"count"] <- cn
		} else { 
			blk[1:(nb-1), "count"] <- block.size
			blk[nb, "count"] <- cn -(nb-1)*block.size
			for(i in 2:nb) blk[i,"start"] <- blk[i-1,"start"] + blk[i-1,"count"]
		}
 
	
		###################  Loop over blocks within chromosome ############################
		
		for(i in 1:nb){  

			# get genotypic data
			bs <- blk[i,"start"]; bc <- blk[i,"count"]
                        geno <- getGenotype(genoData, snp=c(bs, bc), scan=c(1,-1))
			geno <- geno[,sel]  # select samples
		
			##########################  Loop over SNPs within block  ###############################
			
			for(j in 1:bc){

				index <- bs + j - 1
				int.id <- snpID[index]

				# combine with phenotypic data and get complete cases
                                gtype <- geno[j,]
                                pheno.com <- cbind(pheno, gtype)
                                pheno.com <- pheno.com[complete.cases(pheno.com),]
				
				# get MAF
                                maf <- sum(pheno.com$gtype)/(2*length(pheno.com$gtype))
                                maf <- min(maf, 1-maf)
                                mafx <- NA
                                if(chr == XchromCode(genoData)) {  # for X-linked loci
                                        gm <- pheno.com$gtype[pheno.com$sex=="M"]
                                        gf <- pheno.com$gtype[pheno.com$sex=="F"]
                                        mafx <- (sum(gf) + sum(gm)/2)	/ (2*length(gf) + length(gm))	
                                        mafx <- min(mafx, 1-mafx)	
                                }
					
				# decide whether to fit the cox model and do so if criterion is met
                                fit <- !is.na(maf) & maf>0
                                ne <- sum(pheno.com[,event])	
                                if(maf.filter==TRUE & fit==TRUE) {  # additional maf criterion
                                        fit <- maf*(1-maf) > 75/(2*ne)  # use maf rather than mafx for x chromosome loci here
                                }
				#if(chr==XchromCode(genoData)) maf <- mafx  # but record the correct maf in results table	

				if(!fit) {  
					res[k,] <- list(index, int.id, chr, maf, mafx, NA, NA, NA, NA, NA, ne)
					if(!is.null(GxE)){
                                                res2[k,] <- list(index, int.id, chr, maf, mafx, NA, ne, NA, NA)
                                        }
				} else {
					# fit the main effects model
                                        surv <- Surv(pheno.com[,time.to.event], pheno.com[,event])
                                        cph <- tryCatch(coxph( model, data=pheno.com), warning=function(w) TRUE, error=function(e) TRUE)
                                        if(is.logical(cph)) {  # i.e. if an error or warning was issued
                                                res[k,] <- list(index, int.id, chr, maf, mafx, NA, NA, NA, NA, cph, ne)
					} else {
						cphr <- summary(cph)$coef["gtype",][-2]
                                                res[k,] <- list(index, int.id, chr, maf, mafx, cphr["coef"], cphr["se(coef)"], cphr["z"], cphr["Pr(>|z|)"], FALSE, ne)
                                                loglik <- cph$loglik[2]
					}
					# fit the interaction model
					if(!is.null(GxE)){
						cph <- tryCatch(coxph( model2, data=pheno.com), warning=function(w) TRUE, error=function(e) TRUE)
						if(is.logical(cph)) {  # i.e. if an error or warning was issued
							res2[k,] <- list(index, int.id, chr, maf, mafx, cph, ne, NA, NA)
						} else {
							loglik2 <- cph$loglik[2]	
							lrtest <- -2*(loglik-loglik2)
							pval <- 1-pchisq(lrtest,1)
							res2[k,] <- list(index, int.id, chr, maf, mafx, FALSE, ne, lrtest, pval)
						}
					} # end GxE if
						
				}
				#message(paste("snp", k, "of", nsnp))
				k <- k + 1		
										
			}
			#########################  End of loop over SNPs within block  ##########################
			if (verbose) message(paste("block", i, "of", nb, "for chromosome", chr, "and", k-1, "of", nsnp, "SNPs"))
		}	
		######################### End of loop over blocks within chromosome  ####################
	}
	##########################  End of loop over chromosomes  #########################			


    # save the results
    if (!is.null(outfile)) {
      if (allchrom) {
        fileOut <- paste(outfile, "RData", sep=".");
      } else {
        fileOut <- paste(outfile, "chr", paste(chromosome.set[1],chromosome.set[length(chromosome.set)],sep="_"), "RData", sep=".");
      }
      save(res, file=fileOut, compress=TRUE);

      if (!is.null(GxE)) {
        if (allchrom) {
          fileOut <- paste(outfile, "GxE", "RData", sep=".");
        } else {
          fileOut <- paste(outfile, "GxE", "chr", paste(chromosome.set[1],chromosome.set[length(chromosome.set)],sep="_"), "RData", sep=".");
        }
        save(res2, file=fileOut, compress=TRUE);
      }
      
      # save the warnings
      warn <- warnings();
      if (!is.null(warn)) {
        if (allchrom) {
          warnfileOut <- paste(outfile, "warnings", "RData", sep=".");
        } else {
          warnfileOut <- paste(outfile, "chr", paste(chromosome.set[1],chromosome.set[length(chromosome.set)],sep="_"), "warnings", "RData", sep=".");
        }
        save(warn, file=warnfileOut);
      }
      return(invisible(NULL))
    } else {
      if (!is.null(GxE)) {
        return(list("main"=res, "GxE"=res2))
      } else {
        return(res)
      }
    }
}

	

