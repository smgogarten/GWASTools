# Methods for MatrixGenotypeReader

MatrixGenotypeReader <- function(genotype=genotype, snpID=snpID, chromosome=chromosome, position=position, scanID=scanID, ...) {
  new("MatrixGenotypeReader", genotype=genotype, snpID=snpID, chromosome=chromosome, position=position, scanID=scanID, ...)
}

setValidity("MatrixGenotypeReader",
    function (object) {
      # check that dimensions and variables are as expected
      if (length(object@snpID) != nsnp(object)) {
        return("snpID has incorrect dimension")
      }
      if (length(object@chromosome) != nsnp(object)) {
        return("chromosome has incorrect dimension")
      }
      if (length(object@position) != nsnp(object)) {
        return("position has incorrect dimension")
      }      
      if (length(object@scanID) != nscan(object)) {
        return("scanID has incorrect dimension")
      }
      if (!allequal(dim(object@genotype), c(nsnp(object), nscan(object)))) {
        return("genotype has incorrect dimensions")
      }
      
      TRUE
    })

# accessor methods
# index is logical or integer vector of indices to return
setMethod("getSnpID",
          signature(object="MatrixGenotypeReader"),
          function(object, index) {
            var <- object@snpID
            if (missing(index)) var else var[index]
          })

# char=TRUE to return character code
setMethod("getChromosome",
          signature(object="MatrixGenotypeReader"),
          function(object, index, char=FALSE) {
            var <- object@chromosome
            if (!missing(index)) var <- var[index]

            # convert to characters
            if (char) {
              # default is unknown code
              chromChar <- rep("U", length(var))
              autosome <- var %in% object@autosomeCode
              chromChar[autosome] <- as.character(var[autosome])
              xchrom <- var == object@XchromCode & !is.na(var)
              chromChar[xchrom] <- "X"
              ychrom <- var == object@YchromCode & !is.na(var)
              chromChar[ychrom] <- "Y"
              xychrom <- var == object@XYchromCode & !is.na(var)
              chromChar[xychrom] <- "XY"
              mchrom <- var == object@MchromCode & !is.na(var)
              chromChar[mchrom] <- "M"
              var <- chromChar
            }
            var
          })
                        
setMethod("getPosition",
          signature(object="MatrixGenotypeReader"),
          function(object, index) {
            var <- object@position
            if (missing(index)) var else var[index]
          })
                        
setMethod("getScanID",
          signature(object="MatrixGenotypeReader"),
          function(object, index) {
            var <- object@scanID
            if (missing(index)) var else var[index]
          })
 
setMethod("getGenotype",
          signature(object="MatrixGenotypeReader"),
          function(object, snp=snp, scan=scan, ...) {
            if (!missing(snp) & !missing(scan)) {
              # get start and count from snp
              snpstart <- snp[1]
              snpcount <- snp[2]
              if (snpcount == -1) {
                snpend <- nsnp(object)
              } else {
                snpend <- snpstart + snpcount - 1
              }
            
              # get start and count from scan
              scanstart <- scan[1]
              scancount <- scan[2]
              if (scancount == -1) {
                scanend <- nscan(object)
              } else {
                scanend <- scanstart + scancount - 1
              }

              object@genotype[snpstart:snpend, scanstart:scanend]
            } else {
              object@genotype
            }
          })

setMethod("nsnp", "MatrixGenotypeReader",
          function(object) {
            nrow(object@genotype)
          })

setMethod("nscan", "MatrixGenotypeReader",
          function(object) {
            ncol(object@genotype)
          })

setMethod("autosomeCode", "MatrixGenotypeReader",
          function(object) {
            object@autosomeCode
          })
        
setMethod("XchromCode", "MatrixGenotypeReader",
          function(object) {
            object@XchromCode
          })
          
setMethod("YchromCode", "MatrixGenotypeReader",
          function(object) {
            object@YchromCode
          })
              
setMethod("XYchromCode", "MatrixGenotypeReader",
          function(object) {
            object@XYchromCode
          })
              
setMethod("MchromCode", "MatrixGenotypeReader",
          function(object) {
            object@MchromCode
          })
          
setMethod("show", 
          signature(object="MatrixGenotypeReader"),
          function(object) {            
            cat("An object of class", class(object), "\n")
            cat(paste("with", nscan(object), "scans and",
                      nsnp(object), "snps\n"))
          })
