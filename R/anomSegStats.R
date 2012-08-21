#####
# Function to calculate LRR and BAF stats for anomalous segments
#####

anomSegStats <- function(intenData, genoData,
	snp.ids,	# vector of "int.id"s for 'eligible' snps (usually missing.n1<1 and excluding HLA and XTR regions)
	anom,		# data.frame with "scanID", "chromosome", "left.index","right.index","sex"
	centromere,		# data.frame with centromere positions, with variables "chromosome", "left.base", "right.base"
	lrr.cut= -2,   # the function counts the number of eligible LRR values less than this lrr.cut
        verbose=TRUE)
{

	# modified from anom.seg.stats.baf.lrr.v3.r
	# v2 adds calculation of the length of each anomaly
	# v3 adds several new stats and renames some of the old ones; rename some of the input variables in anom
	#     use 'eligible' instead of separate snp.exclude baf and lrr;
	# 	calculate 'nmark.baf' and 'nmark.lrr' for each anomaly regardless of 'method'
	#	allow for missing value codes other than -1 for genotypes
	# v4 adds stats to distiguish wide splits from heterozygous deletions using LOH detection
	#	for each snp, calculate y = min of absolute value of BAF deviations from 0 and 1
	#	then yh = y>0 and genotype call = 0 or 2 (i.e. homozygous)
	#	and ya = y>0 and y<0.25 and genotype call = 0, 2 or missing
	#	for the anom and non-anom regions, get mean of ya and yh and mean of sqrt(ya) and sqrt(yh)
	#	mean of yh called "anom.homo.dev" and "non.anom.homo.dev"
	#	mean of ya called "anom.homo.na.dev" and "non.anom.homo.na.dev"
	#	mean of sqrt(yh) called "anom.homo.sqrt" and "non.anom.homo.sqrt"
	#	mean of sqrt(ya) called "anom.homo.na.sqrt" and non.anom.homo.na.sqrt"
	# v5 make calculation of homo stats optional
	# v6 add MAD for non-anom LRR
	# v7 add stats for detecting terminal anomalies
	# v8 more on terminal anoms
	# v9 add stats for analysis anoms that span the centromere - number of elig markers on each side of centromere and distance they span
	# v10 stats for Cooper metric
	# v11 add check for chromosome matching in anom and cent
	# v12 remove the restriction that any variables calculated by anom stats be excluded from the input data set (they will just be replaced anyway)
	# v13 allows centromere stats for chr 24 (pseudo-autosomal) 
	# v14 require "method" and "anom.id", variables in anom
	# v15 adds a count of the number of eligible LRR markers with value < lrr.cut
	# v16 adds SD and MAD for BAF and LRR for eligible markers within the anomaly
	# v17 adds mean for BAF for eligible markers within the anomaly

	# NOTE: WHEN REVISING THIS FUNCTION TO CREATE A NEW VARIABLE IN THE ANOM DATA.FRAME, DEFINE IT AS NA IN THE SECTION BELOW ("results will be added...")

  # check that intenData has BAF/LRR
  if (!hasBAlleleFreq(intenData)) stop("BAlleleFreq not found in intenData")
  if (!hasLogRRatio(intenData)) stop("LogRRatio not found in intenData")
  
  # check that dimensions of intenData and genoData are equal
  intenSnpID <- getSnpID(intenData)
  genoSnpID <- getSnpID(genoData)
  if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
  intenScanID <- getScanID(intenData)
  genoScanID <- getScanID(genoData)
  if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")
  
	intid <- intenSnpID
	sid <- intenScanID
	indices <- 1:length(intid)
	chrom <- getChromosome(intenData)
	pos <- getPosition(intenData)
	nsnp <- nsnp(intenData)

	# check and supplement centromere
	if(class(centromere)!="data.frame") stop("centromere should be a data.frame")  
	if(!all(is.element(c("chrom","left.base","right.base"), names(centromere)))) stop("names of centromere must include chrom, left.base, right.base")
        centromere$chrom[is.element(centromere$chrom, "X")] <- XchromCode(intenData)
        centromere$chrom[is.element(centromere$chrom, "Y")] <- YchromCode(intenData)
        centromere$chrom[is.element(centromere$chrom, "XY")] <- XYchromCode(intenData)
        centromere$chrom <- as.integer(centromere$chrom)
	centromere <- centromere[, c("chrom","left.base","right.base")]  
	if(!all(is.element(anom$chromosome, centromere$chrom))) stop("centromere must include all 'chromosome' values in anom")
	centromere <- centromere[order(centromere$chrom),]
	if(sum(duplicated(centromere$chrom))!=0) stop("chromosome in centromere must not be duplicated")
	if(!all(is.element(anom$chromosome, centromere$chrom))) stop("chromosome  in centromere does not include all those in anom")
	# get indices of markers adjacent to centromere
	centromere$left.index <- NA
	centromere$right.index <- NA
	for(i in 1:nrow(centromere)){
		left <- centromere$left.base[i] - pos
		right <- pos - centromere$right.base[i]
		chr.sel <- chrom==centromere$chrom[i]
		left.sel <- chr.sel & left>0
		right.sel <- chr.sel & right>0
		if(sum(left.sel)>0) {			# this is needed for the acrocentric chroms with no markers left of the centromere
			 lmin <- min(left[chr.sel & left>0])
			 lind <- indices[chr.sel & left==lmin]
			 centromere$left.index[i] <- lind[length(lind)]	
		}
		if(sum(right.sel)>0) {			# for symmetry even though all chroms should have markers right of the centromere
			 rmin <- min(right[chr.sel & right>0])
			 rind <- indices[chr.sel & right==rmin]
			 centromere$right.index[i] <- rind[length(rind)]	
		}
	}
  
        setw <- which(!is.na(centromere$left.index) & !is.na(centromere$right.index))
        for(i in setw){
            tmpind <- indices[indices>centromere$left.index[i] & indices<centromere$right.index[i]]
            tmpid <- intid[tmpind]
            if(any(is.element(tmpid,snp.ids))) stop("eligible snps in centromere gap")
        }
  
	# check anom
	if(class(anom)!="data.frame") stop("anom should be a data.frame")
	if(!all(is.element(c("scanID", "chromosome", "left.index","right.index","sex","method","anom.id"),names(anom)))) stop("anom does not have required variable names")
	if(!all(is.element(c(anom$left.index,anom$right.index),indices))) stop("left.index and/or right.index are not within range of snp indices")
        ok.chrom <- c(1:22, XchromCode(intenData), XYchromCode(intenData))
	if(!all(anom$chromosome %in% ok.chrom)) stop("chromosome must be integer code for autosomes, X, or XY")

	# check eligible and get selection vector
        eligible <- snp.ids
	if(!all(is.element(eligible,intid))) stop("eligible values must be snpIDs")
	elig <- is.element(intid, eligible)

	# unique samples with anoms
	anom <- anom[order(anom$scanID),]
	snum <- unique(anom$scanID)	
	n <- length(snum)

	# results will be added on to the anom data.frame
	m <- nrow(anom)
	anom$left.method <- rep(NA,m) #  first letter of method (for BAF-LOH merges)
	anom$right.method <- rep(NA,m) # last letter of method
	anom$nmark.all <- rep(NA,m)	#  total number of probes on the array from left.index to right.index inclusive
	anom$nmark.elig <- rep(NA,m)  #  total number of eligible probes on the array from left.index to right.index, inclusive
	anom$left.base <- rep(NA,m)   #  base position corresponding to left.index
	anom$right.base <- rep(NA,m)  #  base position corresponding to right.index
	anom$nbase <- rep(NA,m)	      #  number of bases from left.index to right.index, inclusive
	anom$non.anom.baf.med <- rep(NA, m) # BAF median of non-anomalous segments on all autosomes, using (eligible & is.element(geno,c(1,-1))
	anom$non.anom.lrr.med <- rep(NA, m) # LRR median of non-anomalous segments on all autosomes, using eligible snps
	anom$non.anom.lrr.mad <- rep(NA,m)  # MAD for LRR of non-anomalous segments on all autosomes
	anom$anom.baf.dev.med <- rep(NA,m) # BAF median of deviations of points used to detect anomaly(eligible & is.element(geno,c(1,-1)) from non.anom.baf.med
	anom$anom.baf.dev.5 <- rep(NA,m)  # BAF median of deviations of points used to detect anomaly(eligible & is.element(geno,c(1,-1)) from 0.5
	anom$anom.baf.dev.mean <- rep(NA,m)  # BAF mean of deviations of points used to detect anomaly(eligible & is.element(geno,c(1,-1)) from non.anom.baf.med
	anom$anom.baf.sd <- rep(NA,m)  # SD of BAF deviations from non.anom.baf.med for points used to detect anomaly
	anom$anom.baf.mad <- rep(NA,m) # MAD of BAF deviations from non.anom.baf.med for points used to detect anomaly
	anom$anom.lrr.med <- rep(NA,m)  # LRR median of eligible markers in the anomaly
	anom$anom.lrr.sd <- rep(NA,m)  # SD of LRR for eligible points within the anomaly
	anom$anom.lrr.mad <- rep(NA,m) # MAD of LRR for eligible points within the anomaly
	anom$nmark.baf <- rep(NA,m)  #  number of markers eligible for baf detection
	anom$nmark.lrr <- rep(NA,m)  # number of markers eligible for lrr detection
	# for detecting terminal anomalies
	anom$cent.rel <- rep(NA,m)  # position relative to centromere - left, right, span
	anom$left.most <- rep(NA,m)  # T/F is this the left-most anomaly for this sample-chrom
	anom$right.most <- rep(NA,m)  # T/F is this the right-most anomaly for this sample-chrom
	anom$left.last.elig <- rep(NA,m)  # T/F if it contains the last eligible marker going to the left
	anom$right.last.elig <- rep(NA,m)  # T/F if it contains the last eligible marker going to the right
	anom$left.term.lrr.med <- rep(NA,m)  # median of LRR for all markers from leftmost eligible marker to the left telomere (only calculated for the most distal anom)
	anom$right.term.lrr.med <- rep(NA,m)  # median of LRR for all markers from rightmost eligible marker to the right telomere
	anom$left.term.lrr.n <- rep(NA,m) # sample size for calculating median LRR
	anom$right.term.lrr.n <- rep(NA,m)
	# for analyzing centromere-spanning anomalies
	anom$cent.span.left.elig.n <- rep(NA,m)  # number of eligible markers on the left side of centromere-spanning anomalies
	anom$cent.span.right.elig.n <- rep(NA,m) # number of eligible markers on the right side of centromere-spanning anomalies
	anom$cent.span.left.bases <- rep(NA,m)  # length of anomaly covered by eligible markers on the left side of the centromere
	anom$cent.span.right.bases <- rep(NA,m)  # length of anomaly covered by eligible markers on the right side of the centromere
	anom$cent.span.left.index <- rep(NA,m)  # index of elig marker left-adjacent to centromere
	anom$cent.span.right.index <- rep(NA,m)  # index of elig marker right-adjacent to centromere
	# Cooper metric = sqrt(min(BAF,1-BAF,abs(BAF-M))) for eligible markers, where M is the median BAF over non-anomalous autosomal segments
	anom$bafmetric.anom.mean <- rep(NA,m)
	anom$bafmetric.non.anom.mean <- rep(NA,m)
	anom$bafmetric.non.anom.sd <- rep(NA,m)
	# number of eligible lrr markers with value less than lrr.cut
	anom$nmark.lrr.low <- rep(NA,m)

	# check method values 
	anom$left.method <- substr(anom$method,1,1)
	tmp<- nchar(anom$method)
	anom$right.method <- substr(anom$method, tmp, tmp)
	if(!all(is.element(anom$left.method,c("B","L")) & is.element(anom$right.method,c("F","H","B","L")))) stop("method values should begin with B or L and end with F,H,B,L")

	
	r <- 1  # counter for row number in anom

	for(i in 1:n) {  # for each sample  1:n
		
		# get all data for the ith sample
			dat <- anom[is.element(anom$scanID, snum[i]),]
			sind <- which(is.element(sid, snum[i]))
			geno <- getGenotype(genoData, snp=c(1,-1), scan=c(sind,1))
                        baf <- getBAlleleFreq(intenData, snp=c(1,-1), scan=c(sind,1))
                        lrr <- getLogRRatio(intenData, snp=c(1,-1), scan=c(sind,1))
		# make missing genotype calls equal -1 (regardless of missing value code)
			geno[is.na(geno)] <- -1
		# get vector to select eligible markers not in anomalies 
			k <- nrow(dat)
			ai <- NULL	#indices included in anomalies
			for(j in 1:k) ai <- c(ai, dat$left.index[j]:dat$right.index[j])
			non.baf <- is.element(chrom, 1:22) & is.element(intid, eligible) & !is.element(indices, ai) & is.element(geno, c(1,-1))
			non.lrr <- is.element(chrom, 1:22) & is.element(intid, eligible) & !is.element(indices, ai)
		# calculate median of baf and lrr in non-anomalous regions
			baf.non.med <- median(baf[non.baf], na.rm=TRUE)
			lrr.non.med <- median(lrr[non.lrr], na.rm=TRUE)
		# calculate MAD of lrr in non-anomalous regions
			lrr.non.mad <- mad(lrr[non.lrr], na.rm=TRUE)
		# calculate stats in anomalous segments
		# eligible points used for anomaly detection
			ebaf <- is.element(intid, eligible) & is.element(geno, c(1,-1))
			elrr <- is.element(intid, eligible)
		# Cooper metric
			coop <- sqrt(pmin( (baf-0), (1-baf), abs(baf-baf.non.med) ))
			coop.non.mean <- mean(coop[non.baf], na.rm=TRUE)
			coop.non.sd <- sd(coop[non.baf], na.rm=TRUE)
	
		# loop through anomalies within sample
			for(j in 1:k){ # for each anomaly within a sample

				# MARKER NUMBERS
					anom$nmark.all[r] <- dat$right.index[j] - dat$left.index[j] + 1
					anom$nmark.elig[r] <- sum(is.element(intid,eligible)[dat$left.index[j]:dat$right.index[j]])

				# BAF
					baf.sel <- ebaf & indices>=dat$left.index[j] & indices<=dat$right.index[j] & !is.na(baf)
					baf.dev <- abs(baf[baf.sel]-baf.non.med)
					anom$anom.baf.dev.med[r] <- median(baf.dev)
					anom$anom.baf.dev.mean[r] <- mean(baf.dev)
					baf.dev <- abs(baf[baf.sel]- 0.5)
					anom$anom.baf.dev.5[r] <- median(baf.dev)
					anom$non.anom.baf.med[r] <- baf.non.med
					anom$nmark.baf[r] <- sum(baf.sel)
					anom$anom.baf.sd[r] <- sd(baf.dev)
					anom$anom.baf.mad[r] <- mad(baf.dev)

				# COOPER METRIC
					anom$bafmetric.anom.mean[r] <- mean(coop[baf.sel])
					anom$bafmetric.non.anom.mean[r] <- coop.non.mean
					anom$bafmetric.non.anom.sd[r] <- coop.non.sd

				# LRR
					lrr.sel <- elrr & indices>=dat$left.index[j] & indices<=dat$right.index[j] & !is.na(lrr)
					anom$anom.lrr.med[r] <- median(lrr[lrr.sel])
					anom$non.anom.lrr.med[r] <- lrr.non.med
					anom$non.anom.lrr.mad[r] <- lrr.non.mad
					anom$nmark.lrr[r] <- sum(lrr.sel)
					anom$nmark.lrr.low[r] <- sum(lrr[lrr.sel] < lrr.cut)
					anom$anom.lrr.sd[r] <- sd(lrr[lrr.sel])
					anom$anom.lrr.mad[r] <- mad(lrr[lrr.sel])

				# BASES
					anom$left.base[r] <- pos[dat$left.index[j]]
					anom$right.base[r] <- pos[dat$right.index[j]]
					anom$nbase[r] <- anom$right.base[r] - anom$left.base[r] + 1

				# CENTROMERE
				# find relation to centromere
					chrj <- anom$chromosome[r]
					left <- anom$right.base[r] <= centromere$left.base[centromere$chrom==chrj]
					right <- anom$left.base[r] >= centromere$right.base[centromere$chrom==chrj]
					span <- anom$left.base[r] <= centromere$left.base[centromere$chrom==chrj] & 
						 anom$right.base[r] >= centromere$right.base[centromere$chrom==chrj]
					if(left) anom$cent.rel[r] <- "left"
					if(right) anom$cent.rel[r] <- "right"
					if(span) anom$cent.rel[r] <- "span"
					if(is.na(anom$cent.rel[r])) { stop(paste("i=",i,"j=",j,"\ncent.rel is missing")) }
				# for centromere-spanning, find number of eligible markers (method dependent) and length spanned by elig markers on each side
					if(anom$cent.rel[r]=="span" & !is.element(anom$method[r], c("BAF","LOH")) ) {
						warning(paste("sample",anom$scanID[r],"chromosome", anom$chromosome[r],"spans centromere and has non-standard method", anom$method[r],
						"so centromere check will not be done"))
					}
					if(anom$cent.rel[r]=="span" & is.element(anom$method[r],c("BAF","LOH"))){
						if(is.element(anom$method[r],"LOH")) { elog <- elrr } else { elog <- ebaf }
						left.elig <- which(elog & indices >= anom$left.index[r] & indices <= centromere$left.index[centromere$chrom==chrj])
						right.elig <- which(elog & indices >= centromere$right.index[centromere$chrom==chrj] & indices <= anom$right.index[r])
						nle <- length(left.elig)
						nre <- length(right.elig)
						anom$cent.span.left.elig.n[r] <- nle		# number of eligible markers in anom left of the centromere
						anom$cent.span.right.elig.n[r] <- nre
						anom$cent.span.left.bases[r] <- pos[left.elig[nle]] - pos[left.elig[1]]+ 1   # number of bases covered by eligible markers in the anom left of the centromere  *
						anom$cent.span.right.bases[r] <- pos[right.elig[nre]] - pos[right.elig[1]] + 1
						anom$cent.span.left.index[r] <- left.elig[nle]  # index of eligible marker just left of the centromere (for splitting if necessary)
						anom$cent.span.right.index[r] <- right.elig[1]
					}
						
				# TERMINAL
				# is it the most distal anomaly in this sample-chrom?
					datc <- dat[dat$chromosome==anom$chromosome[r],c("left.index","right.index")]
					anom$left.most[r] <- anom$left.index[r]==min(datc$left.index)
					anom$right.most[r] <- anom$right.index[r]==max(datc$right.index)
				# check whether all markers distal to the breakpoint are ineligible, which depends on method
				# method will be BAF, LOH (for original anoms) or BL, LB for merged anoms (corresponding to BAF left end and LOH right end or vice versa)
					sel.chrom <- chrom==anom$chromosome[r]
					if(anom$left.most[r]){
						anom$left.last.elig[r] <- FALSE
						if(is.element(anom$left.method[r],"L")){ 
							chk.left <- all(!elrr[ sel.chrom & indices < anom$left.index[r] ]) # all distal markers not elig
						}
						if(is.element(anom$left.method[r],"B")){ 
							chk.left <- all(!ebaf[ sel.chrom & indices < anom$left.index[r] ]) # all distal markers not elig
						}
						if(chk.left){	# get median lrr for left-distal markers
							anom$left.last.elig[r] <- TRUE
							chr.ind <- indices[sel.chrom]
							left <- chr.ind[1]; right <- anom$left.index[r]-1
							anom$left.term.lrr.n[r] <- right-left+1
							if(right>=left) anom$left.term.lrr.med[r] <- median(lrr[left:right], na.rm=TRUE) 
						}
					}
					if(anom$right.most[r]){
						anom$right.last.elig[r] <- FALSE
						if(is.element(anom$left.method[r],c("L","H"))){ 
							chk.right <- all(!elrr[ sel.chrom & indices > anom$right.index[r] ]) # all distal markers not elig
						}
						if(is.element(anom$left.method[r],c("B","F"))){ 
							chk.right <- all(!ebaf[ sel.chrom & indices > anom$right.index[r] ]) # all distal markers not elig
						}
						if(chk.right){	# get median lrr for right-distal markers
							anom$right.last.elig[r] <- TRUE
							chr.ind <- indices[sel.chrom]; ni <- length(chr.ind)
							left <- anom$right.index[r]+1; right <- chr.ind[ni]
							anom$right.term.lrr.n[r] <- right-left+1
							if(right>=left) anom$right.term.lrr.med[r] <- median(lrr[left:right], na.rm=TRUE) 
						}
					}

				# increment counter
					#message(paste("r=",r))
					r <- r+1
			}  # end loop over anomalies within a sample

			# print status
			if(verbose & i %% 100 == 0) message(paste("sample", i, "of", n))

		}  # end loop over samples

        # remove "left.method" and "right.method"
        anom$left.method <- NULL
        anom$right.method <- NULL
	return(anom)
}




