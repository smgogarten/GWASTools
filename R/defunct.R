gdsCheckImputedDosage <- function(genoData, snpAnnot, scanAnnot, 
                                  input.files, chromosome,
                                  input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                                  input.dosage=FALSE, block.size=5000,
                                  verbose=TRUE, 
                                  snp.exclude=NULL,
                                  snp.id.start=1,
                                  tolerance=1e-4,
                                  na.logfile=NULL) {

  .Defunct("checkImputedDosageFile")
}

gdsImputedDosage <- function(input.files, gds.filename, chromosome,
                             input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                             input.dosage=FALSE, block.size=5000,
                             snp.annot.filename="dosage.snp.RData",
                             scan.annot.filename="dosage.scan.RData",
                             verbose=TRUE, zipflag="ZIP.max", genotypeDim="snp,scan",
                             scan.df=NULL,
                             snp.exclude=NULL,
                             snp.id.start=1) {

  .Defunct("imputedDosageFile")
}

gdsSetMissingGenotypes <- function(
	parent.file,	# name of the parent GDS file
	new.file,		# name of the new GDS file to create
        regions, # data frame with regions
	sample.include=NULL,	# vector of sampleIDs for samples to include in new.file
        zipflag = "ZIP.max",
        verbose=TRUE) {

    .Defunct("setMissingGenotypes")
}

ncdfAddData <- function(path=".", 
                         ncdf.filename, 
                         snp.annotation,
                         scan.annotation, 
                         sep.type, 
                         skip.num, 
                         col.total, 
                         col.nums, 
                         scan.name.in.file,
                         scan.start.index = 1,
                         diagnostics.filename = "ncdfAddData.diagnostics.RData",
                         verbose = TRUE) {

    .Defunct("createDataFile")
}

ncdfAddIntensity <- function(path=".",
                           ncdf.filename,
                           snp.annotation,
                           scan.annotation, 
                           scan.start.index = 1, 
                           n.consecutive.scans = -1,  
                           diagnostics.filename = "ncdfAddIntensity.diagnostics.RData",
                           verbose = TRUE) {

    .Defunct("createAffyIntensityFile")
}

ncdfCheckGenotype <- function(path=".",
                            ncdf.filename, 
                            snp.annotation, 
                            scan.annotation, 
                            sep.type,
                            skip.num,
                            col.total,
                            col.nums,
                            scan.name.in.file, 
                            check.scan.index,
                            n.scans.loaded,
                            diagnostics.filename = "ncdfCheckGenotype.diagnostics.RData",
                            verbose = TRUE) {
		
    .Defunct("checkGenotypeFile")
}

ncdfCheckIntensity <- function(path=".",
                             intenpath=".",
                             ncdf.filename, 
                             snp.annotation, 
                             scan.annotation, 
                             sep.type,
                             skip.num,
                             col.total,
                             col.nums,
                             scan.name.in.file, 
                             check.scan.index,
                             n.scans.loaded,
                             affy.inten = FALSE,
                             diagnostics.filename = "ncdfCheckIntensity.diagnostics.RData",
                             verbose = TRUE) {

    .Defunct("checkIntensityFile")
}

ncdfCreate <- function(snp.annotation, 
                            ncdf.filename, 
                            variables = "genotype",
                            n.samples = 10, 
                            precision = "double",
                            array.name = NULL,
                            genome.build = NULL) {
    .Defunct("createDataFile")
}

ncdfImputedDosage <- function(input.files, ncdf.filename, chromosome,
                              input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                              input.dosage=FALSE, block.size=5000,
                              snp.annot.filename="dosage.snp.RData",
                              scan.annot.filename="dosage.scan.RData",
                              verbose=TRUE) {
  
  .Defunct("imputedDosageFile")
}

ncdfSetMissingGenotypes <- function(
	parent.file,	# name of the parent netCDF file
	new.file,		# name of the new netCDF file to create
        regions, # data frame with regions
	sample.include=NULL,	# vector of sampleIDs for samples to include in new.file
        verbose=TRUE) {

    .Defunct("setMissingGenotypes")
}
