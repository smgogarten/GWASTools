

ncdfCreate <- function(snp.annotation, 
                            ncdf.filename, 
                            variables = "genotype",
                            n.samples = 10, 
                            precision = "double",
                            array.name = NULL,
                            genome.build = NULL)
{
# v2 changes missing value for R and Theta to -9999
# v3 adds BAlleleFreq and LogRRatio
# v4 changes missing value for varchrom to -1 and for all quant variables to -9999

        # check variables
        stopifnot(all(variables %in% c("genotype", "quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq", "LogRRatio")))
        
        # check that snp annotation has right columns
        stopifnot(all(c("snpID", "chromosome", "position") %in% names(snp.annotation)))

        # make sure all snp annotation columns are integers
        if (!is(snp.annotation$snpID, "integer")) {
          snp.annotation$snpID <- as.integer(snp.annotation$snpID)
          warning(paste("coerced snpID to type integer"))
        }
        if (!is(snp.annotation$chromosome, "integer")) {
          snp.annotation$chromosome <- as.integer(snp.annotation$chromosome)
          warning(paste("coerced chromosome to type integer"))
        }
        if (!is(snp.annotation$position, "integer")) {
          snp.annotation$position <- as.integer(snp.annotation$position)
          warning(paste("coerced position to type integer"))
        }

        # make sure snpID is unique
        stopifnot(length(snp.annotation$snpID) == length(unique(snp.annotation$snpID)))

        # make sure snpID is sorted by chromsome and position
        stopifnot(all(snp.annotation$snpID == sort(snp.annotation$snpID)))
        sorted <- order(snp.annotation$chromosome, snp.annotation$position)
        if (!all(snp.annotation$snpID == snp.annotation$snpID[sorted])) {
          stop("snpID is not sorted by chromosome and position")
        }
      
	# Create the netCDF file and load snp annotation
	# Define dimensions
        snpdim <- dim.def.ncdf("snp","count", snp.annotation$snpID)
        sampledim <- dim.def.ncdf("sample","count",1:n.samples, unlim=TRUE)
        
	# Define variables
        varID<-var.def.ncdf("sampleID","id",dim=sampledim, missval=0, prec="integer")
        varpos <- var.def.ncdf("position","bases",dim=snpdim, missval=-1, prec="integer")
        varchr <- var.def.ncdf("chromosome","id",dim=snpdim, missval=-1, prec="integer")
        vargeno <- var.def.ncdf("genotype", "num_A_alleles", dim=list(snpdim,sampledim), missval=-1, prec="byte")
        varqs <-var.def.ncdf("quality","score",dim=list(snpdim, sampledim), missval=-9999, prec=precision)	#quality score	
        varX <- var.def.ncdf("X", "intensity", dim=list(snpdim,sampledim), missval=-9999, prec=precision)
        varY <- var.def.ncdf("Y", "intensity", dim=list(snpdim, sampledim), missval=-9999, prec=precision)
        varrawX <- var.def.ncdf("rawX", "intensity", dim=list(snpdim, sampledim), missval=-9999, prec=precision)
        varrawY <- var.def.ncdf("rawY", "intensity", dim=list(snpdim, sampledim), missval=-9999, prec=precision)
        varR <- var.def.ncdf("R", "intensity", dim=list(snpdim, sampledim), missval=-9999, prec=precision)
        varTheta <- var.def.ncdf("Theta", "intensity", dim=list(snpdim, sampledim), missval=-9999, prec=precision)
        varBAlleleFreq <- var.def.ncdf("BAlleleFreq", "intensity", dim=list(snpdim, sampledim), missval=-9999, prec=precision)
        varLogRRatio <- var.def.ncdf("LogRRatio", "intensity", dim=list(snpdim, sampledim), missval=-9999, prec=precision)
        
	# Select variables to use
        varlist <- list(varID, varpos, varchr,vargeno,varqs,varX,varY,varrawX,varrawY,varR,varTheta,varBAlleleFreq,varLogRRatio)
        varnames <- c("sampleID", "position", "chromosome","genotype", "quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq", "LogRRatio")
        variables <- c("sampleID", "position", "chromosome", variables)
        varuse <- varlist[is.element(varnames, variables)]
        
	# Create the netCDF file
        genofile <- create.ncdf(ncdf.filename, varuse)
        
	# Add position data
        put.var.ncdf(genofile, varpos, snp.annotation$position)
        put.var.ncdf(genofile, varchr, snp.annotation$chromosome)
        
	# Add attributes
        if (!is.null(array.name)) att.put.ncdf( genofile, 0, "array_name", array.name ) 
        if (!is.null(genome.build)) att.put.ncdf( genofile, 0, "genome_build", genome.build ) 

        close.ncdf(genofile)
        return(invisible(NULL))
}

