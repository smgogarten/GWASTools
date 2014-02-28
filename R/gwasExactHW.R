#When performing genotype counts, this function calls on the HWExact function from the GWASExactHW package by Ian Painter.  
#One thing to note about the HWExact function is that it refers to the alleles as "A" and "a", rather than "A" and "B", so the naming conventions of the two functions don't match exactly, but allele "a" is equivalent to allele "B".


gwasExactHW <- function (genoData,
                         scan.chromosome.filter = NULL,
                         scan.exclude = NULL,
                         geno.counts = TRUE,
                         chromosome.set = NULL,
                         block.size = 5000,                      
                         verbose = TRUE,
                         outfile = NULL) 
{
	

    CountGenotypes <- function(genotypes) {
        nAA <- apply(genotypes, 1, function(x) sum(x == 2, na.rm=TRUE))
        nAa <- apply(genotypes, 1, function(x) sum(x == 1, na.rm=TRUE))
        naa <- apply(genotypes, 1, function(x) sum(x == 0, na.rm=TRUE))
        return(data.frame(nAA, nAa, naa))
    }
    
    
    CountAlleles <- function(gc) {
        nA <- 2 * gc$nAA + gc$nAa
        na <- 2 * gc$naa + gc$nAa
        return(data.frame(nA, na))
    }

    # SNP annotation
    snpID <- getSnpID(genoData);
    position <- getPosition(genoData);
    chrom <- getChromosome(genoData);
    unique_chrom <- unique(chrom);
    nChromosomes <- max(chrom);
    rle_chrom <- rle(chrom);
    rle_chrom2 <- rep(0,nChromosomes);
    rle_chrom2[unique_chrom] <- rle_chrom$lengths;

    # create vector of first SNP number for each chromosome
    gbreaks <- vector(length=nChromosomes+1);
    gbreaks[1] <- 1;
    for(a in 1:nChromosomes)
      gbreaks[a+1] <- sum(rle_chrom2[1:a])+1;

    # chromosome.set check
    if(is.null(chromosome.set)){
      chromosome.set <- unique_chrom;
      allchrom <- TRUE
    } else {
      stopifnot(all(chromosome.set %in% unique_chrom))
      allchrom <- FALSE
    }

    # X chromosome check for sex variable
    if(XchromCode(genoData) %in% chromosome.set & !hasSex(genoData)){
      stop("Sex values for the samples are required to compute MAF for chromosome X SNPs");
    }

    # create a vector of column names for the output
    nv <- c("snpID", "chromosome", "position")
    if(geno.counts) {   
      nv <- append(nv, c("nAA", "nAB", "nBB"))
    }
    nv <- append(nv, c("MAF", "minor.allele", "f"))
    nv <- append(nv, "p.value")

    ############# loop through chromosomes ################
    if(verbose) message("Beginning calculations...")
    res.list <- list() # list to keep results
    for (chr.index in 1:length(chromosome.set)){
      chr <- chromosome.set[chr.index];

      # set blocks of SNPs
      n.chr.snp <- rle_chrom2[chr];
      nblocks <- ceiling(n.chr.snp / block.size);

      # set which samples to keep
      keep <- rep(TRUE, nscan(genoData));

      if(!is.null(scan.chromosome.filter)){
        if(dim(scan.chromosome.filter)[2] != nChromosomes){
          stop("The scan.chromosome.filter matrix must have a column for each chromosome, the column number must match the chromosome number")
        }
        if(!all(row.names(scan.chromosome.filter) == getScanID(genoData))){
          stop("The sample IDs from the sample chromosome filter do not match those from the genoData.")
        }
        keep <- keep & scan.chromosome.filter[,chr];
      }

      if(!is.null(scan.exclude)) 
        keep <- keep & !(getScanID(genoData) %in% scan.exclude);

      if(chr == XchromCode(genoData))
        keep <- keep & (getSex(genoData) == "F");

      if(chr %in% c(YchromCode(genoData), MchromCode(genoData))) {
        keep <- keep & FALSE;
        message("HWE test not valid for Y or M; results will be NA")
      }

      ############## loop through blocks ###################
      for(k in 1:nblocks){
        # get SNP information for the block
        start.snp.block <- (k-1)*block.size  + 1;
        n.snps.block <-  block.size;
        if (start.snp.block + n.snps.block > n.chr.snp)
          n.snps.block <-   n.chr.snp - start.snp.block + 1;
        snp.start.pos <- start.snp.block + gbreaks[chr] - 1;
        snp.end.pos <- snp.start.pos + n.snps.block - 1;

        # get genotypes for the block
        geno <- getGenotype(genoData, snp = c(snp.start.pos, n.snps.block), scan = c(1, -1));
        geno <- geno[,keep];

        # create data.frame of results for this block
        res <- cbind(snpID[snp.start.pos:snp.end.pos], chrom[snp.start.pos:snp.end.pos], position[snp.start.pos:snp.end.pos])
        rownames(res) <- snpID[snp.start.pos:snp.end.pos];

        # count genotypes for this block
        tmpGenotypeCounts <- CountGenotypes(geno);

        # put genotype counts in the results
        if(geno.counts){
          res <- cbind(res, tmpGenotypeCounts);
        }

        # calculate MAF
        alleleCounts <- CountAlleles(tmpGenotypeCounts);
        aFreq <- alleleCounts[,2]/(as.matrix(alleleCounts) %*% cbind(c(1, 1)));
        minor.allele <- ifelse(aFreq > 0.5, "A", "B"); 
        selTmp <- !is.na(aFreq) & (aFreq > 0.5);
        aFreq[selTmp] <- 1 - aFreq[selTmp];

        # calculate inbreeding coefficient
        obs.het <- tmpGenotypeCounts[,2]
        geno.tot <- apply(tmpGenotypeCounts,1,sum)
        exp.het <- 2*aFreq*(1-aFreq)*geno.tot
        f <- 1-(obs.het/exp.het)

        # calculate HW p-val and put in results
        hwePs <- GWASExactHW::HWExact(tmpGenotypeCounts);
        res <- cbind(res,  aFreq, minor.allele, f, hwePs);

        # add results to list
        res.list[[paste(chr.index,k)]] <- res
                
        if(verbose) message(paste("chr",chr," block",k, "of",nblocks, "completed"));
        
      } # block loop
      
    } # chromosome loop
    names(res.list) <- NULL
    res.set <- do.call(rbind, res.list)
    colnames(res.set) = nv;

    # set pvalue to NA for monomorphic SNPs
    res.set$p.value[res.set$MAF == 0] <- NA
    
    # save the results
    if (!is.null(outfile)) {
      if(verbose) message("Saving results...")
      if (allchrom) {
        fileOut <- paste(outfile, "RData", sep=".");
      } else {
        fileOut <- paste(outfile, "chr", paste(chromosome.set[1],chromosome.set[length(chromosome.set)],sep="_"), "RData", sep=".");
      }
      save(res.set, file=fileOut, compress=TRUE);
    
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
      return(res.set)
    }
  }




