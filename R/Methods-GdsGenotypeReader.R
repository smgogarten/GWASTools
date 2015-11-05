# Methods for GdsGenotypeReader

GdsGenotypeReader <- function(filename, genotypeDim, genotypeVar, snpIDvar, scanIDvar, ...) {
  if (missing(filename)) stop("filename is required")
  if (missing(genotypeVar)) genotypeVar <- "genotype"
  if (missing(snpIDvar)) snpIDvar <- "snp.id"
  if (missing(scanIDvar)) scanIDvar <- "sample.id"
  
  # GdsReader does not have ... in its argument
  #tmpobj <- new("GdsReader", GdsReader(filename))
  input.gds <- is(filename, 'gds.class')
  tmpobj <- GdsReader(filename)
  
  # automatic checking for genotypeDim:
  snpDim <- getDimension(tmpobj, snpIDvar)
  scanDim <- getDimension(tmpobj, scanIDvar)
  genoDim <- getDimension(tmpobj, genotypeVar)

  # automatically set genotypeDim if not set
  if (missing(genotypeDim)) {
    if (snpDim == scanDim) {
      genotypeDim <- ""
    } else if (all(genoDim == c(snpDim, scanDim))) {
      genotypeDim <- "snp,scan"
    } else if (all(genoDim == c(scanDim, snpDim))) {
      genotypeDim <- "scan,snp"
    } else {
      genotypeDim <- ""
    }
  }
  
  # karl: in which case ? Which validity method ?
  # problem: if filename is a gds, it is wrong to close it
  # close(tmpobj) # in case it fails the validity method, close and reopen
  
  #tryCatch(new("GdsReader", filename=filename, handler=handler),
  #         error=function(e) if (!input.gds) closefn.gds(handler))
  tryCatch(new("GdsGenotypeReader", tmpobj, genotypeDim=genotypeDim, genotypeVar=genotypeVar,
      snpIDvar=snpIDvar, scanIDvar=scanIDvar, ...),
      error=function(e) {if (!input.gds) close(tmpobj)
        stop(e)})
}

setValidity("GdsGenotypeReader",
    function (object) {
      # check that dimensions and variables are as expected
      
      # check that variables are snpID, chromosome, position, scanID, genotype
      if (!hasVariable(object, object@snpIDvar)) {
        return(paste("variable", object@snpIDvar, "not found in", object@filename))
      }
      if (!hasVariable(object, object@chromosomeVar)) {
        return(paste("variable", object@chromosomeVar, "not found in", object@filename))
      } 
      if (!hasVariable(object, object@positionVar)) {
        return(paste("variable", object@positionVar, "not found in", object@filename))
      }
      if (!hasVariable(object, object@scanIDvar)) {
        return(paste("variable", object@scanIDvar, "not found in", object@filename))
      }
      if (!hasVariable(object, object@genotypeVar)) {
        return(paste("variable", object@genotypeVar, "not found in", object@filename))
      }
      
      # check that chromosome and position have same dimension as snpID
      snpDim <- getDimension(object, object@snpIDvar)
      chrDim <- getDimension(object, object@chromosomeVar)
      if (chrDim != snpDim) {
        return(paste("variable", object@chromosomeVar, "has incorrect dimension"))
      }
      posDim <- getDimension(object, object@positionVar)
      if (posDim != snpDim) {
        return(paste("variable", object@positionVar, "has incorrect dimension"))
      }
      
      
      #  check that genotype has dimensions [snpID,scanID] or [scanID,snpID]
      if (!(object@genotypeDim %in% c("snp,scan", "scan,snp")))
        return("genotype order is not specified: 'snp,scan' or 'scan,snp'")
      
      scanDim <- getDimension(object, object@scanIDvar)
      genoDim <- getDimension(object, object@genotypeVar)
      if (length(genoDim) != 2) {
        return(paste("variable", object@genotypeVar, "has incorrect dimension"))
      } else if ((object@genotypeDim == "snp,scan" & !all(genoDim == c(snpDim, scanDim))) | 
               (object@genotypeDim == "scan,snp" & !all(genoDim == c(scanDim, snpDim)))) {
        return(paste("variable", object@genotypeVar, "has incorrect dimensions"))
      }
      TRUE
    })

         
# accessor methods
# index is logical or integer vector of indices to return
.logicalIndex <- function(index, dim) {
    if (is.logical(index)) {
        if (length(index) != dim) stop("index has incorrect dimension")
        index
    } else if (is.numeric(index)) {
        x <- rep(FALSE, dim)
        x[index] <- TRUE
        x
    } else {
        stop("index must be logical or integer")
    }
}

setMethod("getVariable",
          signature(object="GdsGenotypeReader"),
          function(object, varname, index, ...) {
              if (missing(index)) {
                  callNextMethod(object, varname, ...)
              } else {
                  dim <- getDimension(object, varname)
                  if (is.list(index)) {
                      sel <- list()
                      for (i in 1:length(index)) {
                          sel[[i]] <- .logicalIndex(index[[i]], dim[i])
                      }
                  } else {
                      sel <- .logicalIndex(index, dim)
                  }
                  callNextMethod(object, varname, sel=sel, ...)
              }                    
          })


setMethod("getSnpID",
          signature(object="GdsGenotypeReader"),
          function(object, ...) {
            getVariable(object, object@snpIDvar, ...)
          })

# char=TRUE to return character code
setMethod("getChromosome",
          signature(object="GdsGenotypeReader"),
          function(object, index, char=FALSE) {
            if (missing(index))
              var <- getVariable(object, object@chromosomeVar)
            else
              var <- getVariable(object, object@chromosomeVar, index)

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
          signature(object="GdsGenotypeReader"),
          function(object, ...) {
            getVariable(object, object@positionVar, ...)
          })
                         
setMethod("getAlleleA",
          signature(object="GdsGenotypeReader"),
          function(object, ...) {
            var <- getVariable(object, object@alleleVar, ...)
            substr(var, 1, regexpr("/", var, fixed=TRUE)-1)
          })
                           
setMethod("getAlleleB",
          signature(object="GdsGenotypeReader"),
          function(object, ...) {
            var <- getVariable(object, object@alleleVar, ...)
            substr(var, regexpr("/", var, fixed=TRUE)+1, nchar(var))
          })
                        
setMethod("getScanID",
          signature(object="GdsGenotypeReader"),
          function(object, ...) {
            getVariable(object, object@scanIDvar, ...)
          })


.startCountToIndex <- function(start, count, total) {
    if (count == -1) {
        seq(start, total)
    } else {
        seq(start, length.out=count)
    }
}
   
## add names to genotype matrix
.addNamesGds <- function(object, var, snp, scan) {
    if (is.matrix(var)) {
        if (object@genotypeDim == "scan,snp") {
            dimnames(var) <- list(getScanID(object, index=scan), getSnpID(object, index=snp))
        } else if (object@genotypeDim == "snp,scan") {
            dimnames(var) <- list(getSnpID(object, index=snp), getScanID(object, index=scan))
        }
    }
    var
}

## return matrix as snp,scan unless transpose=TRUE
.returnSnpScan <- function(object, var, transpose) {
    ## check about returning transposed results based on @genotyepDim:
    ## scan,snp should by default return the transpose of what's in the gds file
    if (!transpose & is.matrix(var) & object@genotypeDim == "scan,snp") return(t(var))
    ## snp,scan should by default return what's in the gds file
    if (transpose & is.matrix(var) & object@genotypeDim == "snp,scan") return(t(var))
    var
}
 
setMethod("getGenotype",
          signature(object="GdsGenotypeReader"),
          function(object, snp=c(1,-1), scan=c(1,-1), drop=TRUE, use.names=FALSE, transpose=FALSE, ...) {

              ## check if we need to switch snp/scan here for a transposed genotype file
              if (object@genotypeDim == "scan,snp") {
                  start <- c(scan[1], snp[1])
                  count <- c(scan[2], snp[2])
              } else if (object@genotypeDim == "snp,scan"){
                  start <- c(snp[1], scan[1])
                  count <- c(snp[2], scan[2])
              }
              var <- getVariable(object, object@genotypeVar, start=start, count=count, drop=FALSE, ...)
              ## set missing values to NA
              var[var == object@missingValue] <- NA

              if (use.names) {
                  snp.ind <- .startCountToIndex(snp[1], snp[2], nsnp(object))
                  scan.ind <- .startCountToIndex(scan[1], scan[2], nscan(object))
                  var <- .addNamesGds(object, var, snp=snp.ind, scan=scan.ind)
              }

              if (drop) var <- drop(var)
              
              ## return matrix as snp,scan unless transpose=TRUE
              .returnSnpScan(object, var, transpose)
          })

## order matrix (which is a subset) by original selection index
.orderBySelection <- function(object, var, snp, scan) {
    if (is.numeric(snp) & !identical(snp, sort(snp))) {
        snp.ind <- na.omit(match(1:nsnp(object), snp))
        if (object@genotypeDim == "scan,snp") {
            var <- var[,snp.ind]
        } else if (object@genotypeDim == "snp,scan"){
            var <- var[snp.ind,]
        }
    }
    if (is.numeric(scan) & !identical(scan, sort(scan))) {
        scan.ind <- na.omit(match(1:nscan(object), scan))
        if (object@genotypeDim == "scan,snp") {
            var <- var[scan.ind,]
        } else if (object@genotypeDim == "snp,scan"){
            var <- var[,scan.ind]
        }
    }
    var
}

setMethod("getGenotypeSelection",
          signature(object="GdsGenotypeReader"),
          function(object, snp=NULL, scan=NULL, snpID=NULL, scanID=NULL, drop=TRUE, use.names=TRUE,
                   order=c("file", "selection"), transpose=FALSE, ...) {

              order <- match.arg(order)

              if (!is.null(snpID)) {
                  if (!is.null(snp)) stop("cannot specify both snp and snpID")
                  snp <- getSnpID(object) %in% snpID
              }
              if (!is.null(scanID)) {
                  if (!is.null(scan)) stop("cannot specify both scan and scanID")
                  scan <- getScanID(object) %in% scanID
              }
              
              if (is.null(snp)) snp <- rep(TRUE, nsnp(object))
              if (is.null(scan)) scan <- rep(TRUE, nscan(object))
              
              ## check if we need to switch snp/scan here for a transposed genotype file
              if (object@genotypeDim == "scan,snp") {
                  sel <- list(scan, snp)
              } else if (object@genotypeDim == "snp,scan"){
                  sel <- list(snp, scan)
              }
              var <- getVariable(object, object@genotypeVar, index=sel, drop=FALSE, ...)
              ## set missing values to NA
              var[var == object@missingValue] <- NA

              if (use.names) var <- .addNamesGds(object, var, snp, scan)

              if (order == "selection") var <- .orderBySelection(object, var, snp, scan)
              
              if (drop) var <- drop(var)
              
              ## return matrix as snp,scan unless transpose=TRUE
              .returnSnpScan(object, var, transpose)
          })


setMethod("nsnp", "GdsGenotypeReader",
          function(object) {
            getDimension(object, object@snpIDvar)
          })

setMethod("nscan", "GdsGenotypeReader",
          function(object) {
            getDimension(object, object@scanIDvar)
          })

setMethod("autosomeCode", "GdsGenotypeReader",
          function(object) {
            object@autosomeCode
          })
   
setMethod("XchromCode", "GdsGenotypeReader",
          function(object) {
            object@XchromCode
          })
          
setMethod("YchromCode", "GdsGenotypeReader",
          function(object) {
            object@YchromCode
          })
              
setMethod("XYchromCode", "GdsGenotypeReader",
          function(object) {
            object@XYchromCode
          })
              
setMethod("MchromCode", "GdsGenotypeReader",
          function(object) {
            object@MchromCode
          })
          
