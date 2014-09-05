gdsSetMissingGenotypes <- function(
	parent.file,	# name of the parent GDS file
	new.file,		# name of the new GDS file to create
        regions, # data frame with regions
	sample.include=NULL,	# vector of sampleIDs for samples to include in new.file
        zipflag = "ZIP.max",
        verbose=TRUE) {

  stopifnot(all(c("scanID", "chromosome", "left.base", "right.base", "whole.chrom") %in% names(regions)))

  gds.old <- GdsGenotypeReader(parent.file)
  snpID <- getSnpID(gds.old)
  chrom <- getChromosome(gds.old)
  pos <- getPosition(gds.old)
  scanID <- getScanID(gds.old)
  if (is.null(sample.include)) {
    sample.include <- scanID
  } else {
    stopifnot(all(sample.include %in% scanID))
  }

  gfile <- createfn.gds(new.file)

  add.gdsn(gfile, "sample.id", sample.include, compress=zipflag, closezip=TRUE)
  add.gdsn(gfile, "snp.id", snpID, compress=zipflag, closezip=TRUE)
  if (hasVariable(gds.old, "snp.rs.id")) {
    add.gdsn(gfile, "snp.rs.id", getVariable(gds.old, "snp.rs.id"),
             compress=zipflag, closezip=TRUE)
  }
  if (hasVariable(gds.old, "snp.allele")) {
    add.gdsn(gfile, "snp.allele", getVariable(gds.old, "snp.allele"),
             compress=zipflag, closezip=TRUE)
  }
  add.gdsn(gfile, "snp.position", pos, compress=zipflag, closezip=TRUE)
  add.gdsn(gfile, "snp.chromosome", chrom, storage="uint8",
           compress=zipflag, closezip=TRUE)
  put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "autosome.start", min(autosomeCode(gds.old)))
  put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "autosome.end", max(autosomeCode(gds.old)))
  put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "X", XchromCode(gds.old))
  put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "XY", XYchromCode(gds.old))
  put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "Y", YchromCode(gds.old))
  put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "M", MchromCode(gds.old))
  put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "MT", MchromCode(gds.old))
  sync.gds(gfile)

  gGeno <- add.gdsn(gfile, "genotype", storage="bit2",
                    valdim=c(length(snpID), length(sample.include)))
  put.attr.gdsn(gGeno, "snp.order")

  ind.old <- which(scanID %in% sample.include)
  for (ind.new in 1:length(sample.include)) {
    if (verbose & ind.new%%10==0) message(paste("sample", ind.new, "of", length(sample.include)))
    geno <- getGenotype(gds.old, snp=c(1,-1), scan=c(ind.old[ind.new],1))
    reg <- regions[regions$scanID %in% scanID[ind.old[ind.new]],]
    if (nrow(reg) > 0) {
      for (a in 1:nrow(reg)) {
        if (reg$whole.chrom[a]) {
          geno[chrom == reg$chromosome[a]] <- NA
        } else {
          geno[chrom == reg$chromosome[a] & reg$left.base[a] <= pos & pos <= reg$right.base[a]] <- NA
        }
      }
    }
    geno[!(geno %in% c(0,1,2))] <- 3
    write.gdsn(gGeno, geno, start=c(1,ind.new), count=c(-1,1))
  }

  closefn.gds(gfile)
  close(gds.old)
  cleanup.gds(new.file)
  return(invisible(NULL))
}
