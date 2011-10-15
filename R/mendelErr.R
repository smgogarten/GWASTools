########################################################################
#
# Return a "mendelList" object
#
# men.list is a list of lists
# first level list is families
# second level list is offspring within families who have one or both parents genotyped
# within the second level is a data.frame with columns offspring, father, mother (scanID)
# when replicates of the same subject.id occur,
# this data.frame has multiple rows representing all combinations of scanIDs
#
########################################################################
mendelList <- function(familyid,
                            offspring,
                            father,
                            mother,
                            sex,
                            scanID)
{
  # Perform some checks
  if(missing(familyid)){
    stop("familyid is missing - this is a vector of family IDs") }
  if(missing(offspring)){
    stop("offspring is missing - this is a vector of offspring IDs") }
  if(missing(father)){
    stop("father is missing - this is a vector of father IDs") }
  if(missing(mother)){
    stop("mother is missing - this is a vector of mother IDs") }
  if(missing(sex)){
    stop("sex is missing - this is a vector of sexes - M/F") }
  if(missing(scanID)){
    stop("scanID is missing - this is a vector of unique sample identifiers") }
  if(!all(sex %in% c("M","F"))){
    stop("sex values must be M/F") }

  if(length(scanID) != length(familyid)){
    stop("The family ID vector does not match the scanID vector in length") }
  if(length(scanID) != length(offspring)){
    stop("The offspring ID vector does not match the scanID vector in length") }
  if(length(scanID) != length(father)){
    stop("The father ID vector does not match the scanID vector in length") }
  if(length(scanID) != length(mother)){
    stop("The mother ID vector does not match the scanID vector in length") }
  if(length(scanID) != length(sex)){
    stop("The sex vector does not match the scanID vector in length") }
  
  # generate "gender", 1 -- male, 2 -- female
  gender <- rep(1, length(sex))
  gender[sex=="F"] <- 2
   
  rv <- NULL
  # family loop
  for (fid in levels(factor(familyid))){  # this is looping through all unique family IDs
    frv <- NULL
    fflag <- familyid==fid  # flag which samples are part of this family

    # offspring loop
    for (childid in levels(factor(offspring[fflag]))){  # this is looping through all offspring in this family
      crv <- NULL

      # get index(es) for this child in offspring ID vector
      for (childindex in which(offspring==childid & fflag)){
        # get index(es) for the father in the offspring ID vector
        fatset <- which(offspring==father[childindex] & fflag)
        # get index(es) for the mother in the offspring ID vector
        motset <- which(offspring==mother[childindex] & fflag)  

        if ((length(fatset)==0) & (length(motset)>0))
          fatset <- -1
        if ((length(motset)==0) & (length(fatset)>0))
          motset <- -1
        
        for (fatindex in fatset){
          for (motindex in motset){
            # offspring information
            res <- scanID[childindex];
            # father information
            if (fatindex > 0){
              if (gender[fatindex] != 1)
                stop(paste("Subject", fatindex, "is a father but is not male."))
              res <- append(res,scanID[fatindex])
            } else{
              res <- append(res,-1)
            }            
            # mother information
            if (motindex > 0){
              if (gender[motindex] != 2)
                stop(paste("Subject", motindex, "is a mother but is not female."))
              res <- append(res,scanID[motindex])
            } else{
              res <- append(res,-1)
            }

            # add offspring/father/mother trio to the results for this child
            crv <- rbind(crv,data.frame(offspring=res[1], father=res[2], mother=res[3]))
          }
        }
      }
      # add this child's results to the family's
      if (!is.null(crv)){
        frv <- append(frv, list(crv))
        # name is offspring ID
        names(frv)[length(frv)] <- childid
      }
    }	# loop for offspring

    # add this family's results to the rest
    if (!is.null(frv)){
      rv <- append(rv, list(frv))
      # name is the family ID
      names(rv)[length(rv)] <- fid
    }
  }	# loop for family

  if (!is.null(rv)) class(rv) <- "mendelList"
  return(rv)
}



########################################################################
#
# Convert a "mendelList" object into a data frame
#
#  mendel.list -- a "mendelList" object
#
########################################################################
mendelListAsDataFrame <- function(mendel.list)
{
	stopifnot(class(mendel.list)=="mendelList")
	rv <- NULL
        # for all families
	for (famidx in 1:length(mendel.list)){
          # for all children in that family
          for (childidx in 1:length(mendel.list[[famidx]]))
            rv <- rbind(rv, mendel.list[[famidx]][[childidx]])
	}
	return(rv)
}



###############################################################################################
#
# Return a list summarizing Mendelian error,
#	list(trios, all.trios, snp)
#
#  genoData   -- GenotypeData object
#  mendel.list	            -- a "mendelList" object specifying trios
#  snp.exclude	            -- subset of SNPs to exclude, an integer vector
#  error.by.snp	            -- output Mendelian errors per SNP, by default(TRUE)
#  error.by.snp.trio        -- output Mendelian errors per SNP for EACH TRIO, by default(FALSE)
#  verbose		    -- whether or not to show progress, by default(TRUE)
#  outfile                  -- a character string giving a name to save the output as
#
###############################################################################################
mendelErr <- function(genoData,
                       mendel.list,
                       snp.exclude = NULL,
                       error.by.snp = TRUE,
                       error.by.snp.trio = FALSE,
                       verbose = TRUE,
                       outfile = NULL)
{

  # get chromosome numbers
  chr <- getChromosome(genoData, char=TRUE)

  # check
  if (!is(mendel.list, "mendelList"))
    stop("the argument 'mendel.list' must be a mendelList object")
  stopifnot(hasSex(genoData))
  
  # logical vector for excluding SNPs
  if (!is.null(snp.exclude)) {
    snpID <- getSnpID(genoData)
    SNPsubset <- !(snpID %in% snp.exclude)
  }
  
  # for Mendelian errors, check autosome, X, XY and Y
  if (!is.null(snp.exclude)) chr <- chr[SNPsubset]
  chr.auto <- is.element(chr, 1:22)
  chr.x <- is.element(chr, "X")
  chr.xy <- is.element(chr, "XY")
  chr.y <- is.element(chr, "Y")
  chr.m <- is.element(chr, "M")
  chr.men <- chr.auto | chr.x | chr.xy | chr.y

  # divide chr into 5 cases
  # 0 -- autosomes (1..22) and XY
  # 1 -- X
  # 2 -- Y
  # 3 -- mtDNA
  # 4 -- missing
  nchr <- rep(4, length(chr)) # default is missing
  nchr[chr.auto] <- 0
  nchr[chr.xy] <- 0
  nchr[chr.x] <- 1
  nchr[chr.y] <- 2
  nchr[chr.m] <- 3

  #  returns an integer vector:
  #  0 -- unable to judge
  #  1 -- SNP is not correct
  #  2 -- SNP is correct
  # mat, dimen: offspring, father, mother
  MMat <- c(   # 64, autosome and XY
    # missing,  BB,        AB,       AA, father
    0,0,0,0 , 0,2,2,1 , 0,2,2,2 , 0,1,2,2 , # missing, mother, every group of 4 values corresponds to the genotype of offspring (missing, BB, AB, AA)
    0,2,2,1 , 0,2,1,1 , 0,2,2,1 , 0,1,2,1 , # BB
    0,2,2,2 , 0,2,2,1 , 0,2,2,2 , 0,1,2,2 , # AB
    0,1,2,2 , 0,1,2,1 , 0,1,2,2 , 0,1,1,2 ) # AA
  MMatXM <- c(  # chr. X male offspring
    # missing,  BB,       AB,       AA, father
    0,0,0,0 , 0,0,0,0 , 0,0,0,0 , 0,0,0,0 , # missing, mother
    0,2,0,1 , 0,2,0,1 , 0,2,0,1 , 0,2,0,1 , # BB, the case that offspring
    0,2,0,2 , 0,2,0,2 , 0,2,0,2 , 0,2,0,2 , # AB  is AB is not detected
    0,1,0,2 , 0,1,0,2 , 0,1,0,2 , 0,1,0,2 ) # AA  here.
  MMatXF <- c(  # chr. X female offspring
    # missing,  BB,       AB,       AA, father
    0,0,0,0 , 0,2,2,1 , 0,0,0,0 , 0,1,2,2 , # missing, mother
    0,2,2,1 , 0,2,1,1 , 0,2,2,1 , 0,1,2,1 , # BB
    0,2,2,2 , 0,2,2,1 , 0,2,2,2 , 0,1,2,2 , # AB
    0,1,2,2 , 0,1,2,1 , 0,1,2,2 , 0,1,1,2 ) # AA
  MMatYM <- c( # chr. Y male offspring
    # missing,  BB,       AB,       AA, father
    0,0,0,0 , 0,2,1,1 , 0,1,2,1 , 0,1,1,2 , # missing, mother
    0,0,0,0 , 0,2,1,1 , 0,1,2,1 , 0,1,1,2 , # BB
    0,0,0,0 , 0,2,1,1 , 0,1,2,1 , 0,1,1,2 , # AB
    0,0,0,0 , 0,2,1,1 , 0,1,2,1 , 0,1,1,2 ) # AA
  MMatYF <- rep(0, 64)
  MMatmtDNA <- c( # mtDNA
    # missing,  BB,       AB,       AA, father
    0,0,0,0 , 0,0,0,0 , 0,0,0,0 , 0,0,0,0 , # missing, mother
    0,2,1,1 , 0,2,1,1 , 0,2,1,1 , 0,2,1,1 , # BB
    0,1,2,1 , 0,1,2,1 , 0,1,2,1 , 0,1,2,1 , # AB
    0,1,1,2 , 0,1,1,2 , 0,1,1,2 , 0,1,1,2 ) # AA
  MMatChrMiss <- rep(0, 64)
  # therefore, totally 4*4*4*10 = 640 cases
  Mat <- c( MMat, MMatXM, MMatYM, MMatmtDNA, MMatChrMiss,
    MMat, MMatXF, MMatYF, MMatmtDNA, MMatChrMiss )

  # index for checking progress
  if(verbose){
    iInfo <- 1 }

  # get Sample IDs
  sampleID <- getScanID(genoData)

  # generate "gender", 1 -- male, 2 -- female
  sex <- getSex(genoData)
  gender <- rep(1, length(sex))
  gender[sex=="F"] <- 2

  # define result data.frame
  if(verbose) message("Preparing data structures...")
  all.trios <- NULL

  # set a vector of null genotypes for missing father or mother
  miss.geno <- rep(0, nsnp(genoData))
  if (!is.null(snp.exclude)){
    miss.geno <- miss.geno[SNPsubset]
  }

  # result -- # of snps for checking and # of snps with errors
  if(error.by.snp){
    snpID <- getSnpID(genoData)
    if (!is.null(snp.exclude)){
      snpID <- snpID[SNPsubset]
    }
    check.cnt <- rep(0, length(miss.geno)); names(check.cnt) <- snpID
    error.cnt <- rep(0, length(miss.geno)); names(error.cnt) <- snpID
    snp <- list(check.cnt = check.cnt, error.cnt = error.cnt)
  } else{
    snp <- NULL
  }

  
  # Family loop
  for (famidx in 1:length(mendel.list)){
    familyid <- names(mendel.list)[famidx]

    # Child (Subject) loop
    for (childidx in 1:length(mendel.list[[famidx]])){
      childid <- names(mendel.list[[famidx]])[childidx]

      # get the trios for this family-child pair
      trio <- mendel.list[[famidx]][[childidx]]

      if(error.by.snp){
        # indicators for SNPs checked and  errors
        check.ind <- rep(FALSE, length(snp$check.cnt))
        error.ind <- rep(FALSE, length(snp$error.cnt))
      }

      # loop through Samples (multiple if there are duplicates of a subject)
      for (i in 1:dim(trio)[1]){
        # each trio
        rtrio <- trio[i, ]

        # genotype of offspring
        gchild <- getGenotype(genoData, snp=c(1,-1), scan=c(which(sampleID == rtrio$offspring), 1))
        gchild[is.na(gchild)] <- -1
        gchild <- gchild + 1
        if (!is.null(snp.exclude)){
          gchild <- gchild[SNPsubset]
        }

        # genotype of father
        if (rtrio$father > 0){
          gfather <- getGenotype(genoData, snp=c(1,-1), scan=c(which(sampleID == rtrio$father), 1))
          gfather[is.na(gfather)] <- -1
          gfather <- gfather + 1
          if (!is.null(snp.exclude)){
            gfather <- gfather[SNPsubset]
          }
        }else{
          gfather <- miss.geno
        }

        # genotype of mother
        if (rtrio$mother > 0){
          gmother <- getGenotype(genoData, snp=c(1,-1), scan=c(which(sampleID == rtrio$mother), 1))
          gmother[is.na(gmother)] <- -1
          gmother <- gmother + 1
          if (!is.null(snp.exclude)){
            gmother <- gmother[SNPsubset]
          }
        }else{
          gmother <- miss.geno
        }

        # divide chr and gender into 10 cases
        mask <- nchr + (gender[which(sampleID == rtrio$offspring)]-1)*5

        # match
        m <- Mat[gchild + gfather*4 + gmother*16 + mask*64 + 1]

        # updating snp - we only want to count a subject once (don't count each sample for a subject)
        if(error.by.snp){
          check.ind <- check.ind | (m > 0)
          error.ind <- error.ind | (m == 1)
        }

        # updating trio
        # total error counts and summary info
        r <- data.frame(fam.id = familyid, child.id = childid,
                        child.scanID = rtrio$offspring, father.scanID = rtrio$father, mother.scanID = rtrio$mother,
                        Men.err.cnt = sum(m==1 & chr.men), Men.cnt = sum(m>0 & chr.men),
                        mtDNA.err = sum(m==1 & chr.m), mtDNA.cnt = sum(m>0 & chr.m),
                        stringsAsFactors = FALSE)

        # error counts by chromosome
        for (j in c(1:22, "X", "XY", "Y")){
          r[[paste("chr", j, sep="")]] <- sum(m==1 & is.element(chr, j))
        }

        # bind to other results
        all.trios <- rbind(all.trios, r)

        # info
        if (verbose){
          message(date(), "\t[ ", iInfo, " ] family: ", familyid,", child: ", childid, ", ", i)
          iInfo <- iInfo + 1
        }
      } # end samples loop


      # update family-child snp info
      if(error.by.snp){
        snp$check.cnt <- snp$check.cnt + ifelse(check.ind,1,0)
        snp$error.cnt <- snp$error.cnt + ifelse(error.ind,1,0)
        if(error.by.snp.trio){
          snp[[paste(familyid, childid, sep=".")]] <- ifelse(error.ind,1,0)
        }
      }
    } # end child (subject) loop
  } # end family loop


  # trios, average values for duplicate samples
  if(verbose) message("Averaging over duplicate trios...")
  if (!is.null(all.trios)){
    d <- paste(all.trios$fam.id, all.trios$child.id)
    # get unique family-child pairs
    id <- unique(d)
    trios <- all.trios[1:length(id), !(names(all.trios) %in% c("child.scanID", "father.scanID", "mother.scanID"))]
    n <- dim(trios)[2]; Ind <- 1
    for (i in id){
      ti <- all.trios[d==i,]
      trios[Ind, 1:2] <- ti[1, 1:2]
      # average counts over all duplicate samples for that child
      for (j in 3:n){
        trios[Ind, j] <- mean(ti[, j+3], na.rm=T)
      }
      Ind <- Ind + 1
    }
  }else{
    trios <- NULL
  }

  rv <- list("trios"=trios, "all.trios"=all.trios, "snp"=snp)
  class(rv) <- "mendelClass"

  if (!is.null(outfile)) {
    # save the output
    if(verbose) message("Saving output...")
    fileOut <- paste(outfile, "RData", sep=".");
    save(rv, file=fileOut, compress=TRUE);

    # save the warnings
    warn <- warnings();
    if (!is.null(warn)) {
      warnfileOut <- paste(outfile, "warnings", "RData", sep=".");
      save(warn, file=warnfileOut);
    }
    return(invisible(NULL))
  } else {
    return(rv)
  }
}
