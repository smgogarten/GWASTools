\name{plinkToNcdf}
\alias{plinkToNcdf}
\title{Create a netCDF file and annotation suitable for use in GWASTools
  from PLINK files}
\description{
\code{plinkToNcdf} creates a netCDF file and scan and SNP annotation
objects from a set of ped and map files.
}
\usage{
plinkToNcdf(pedFile, mapFile, nSamples,
  ncdfFile, snpAnnotFile, scanAnnotFile,
  ncdfXchromCode=23, ncdfXYchromCode=24, ncdfYchromCode=25,
  ncdfMchromCode=26, ncdfUchromCode=27,
  pedMissingCode=0, verbose=TRUE)
}
\arguments{
  \item{pedFile}{PLINK ped file.}
  \item{mapFile}{PLINK map file.  Columns should be chromosome,
    rsID, map distance (not used, but included in output annotation),
    and base-pair position.  If this is an extended map file (.bim),
    columns 5 and 6 will be used to encode allele A and allele B.}
  \item{nSamples}{Number of samples in the ped file.}
  \item{ncdfFile}{Output netCDF file.}
  \item{snpAnnotFile}{Output .RData file for storing a \code{\link{SnpAnnotationDataFrame}}.}
  \item{scanAnnotFile}{Output .RData file for storing a
    \code{\link{ScanAnnotationDataFrame}}.}
  \item{ncdfXchromCode}{Integer value used to represent the X chromosome
    in the netCDF file.  Values of "X" or "23" in the map file are converted to this code.}
  \item{ncdfXYchromCode}{Integer value used to represent the
    pseudoautosomal region of the X and Y chromsomes in the netCDF
    file.  Values of "XY" or "25" in the map file are converted to this code.}
  \item{ncdfYchromCode}{Integer value used to represent the Y chromosome
    in the netCDF file.  Values of "Y" or "24" in the map file are converted to this code.}
  \item{ncdfMchromCode}{Integer value used to represent mitochondrial SNPs
    in the netCDF file.  Values of "MT" or "26" in the map file are converted to this code.}
  \item{ncdfUchromCode}{Integer value used to represent unknown chromosome
    in the netCDF file.  Any values in the map file not in (1:26, "X",
    "Y", "XY", "MT") are converted to this code.}
  \item{pedMissingCode}{Missing genotype code in the ped file.}
  \item{verbose}{logical for whether to show progress information.}
}
\details{
  The netCDF file stores genotype data in byte format, so the PLINK
  genotype is converted to number of A alleles (0, 1, 2, or missing).
  The definitions of A and B alleles may be provided in the map file (column
  5=A, column 6=B).  Otherwise, A and B definitions will be based on the
  order alleles are encountered in the ped file.  (Note that converting
  between ped/map format and bed/bim/fam format in PLINK will not always
  preserve the order of chromosomes, so use caution when matching a bim
  file to a ped file!)
  
  The first six columns of the ped file will be converted to a
  \code{\link{ScanAnnotationDataFrame}}. 
  If the Individual ID (second column of the ped file) contains unique
  integers, then this column will be used for scanID.  Otherwise, an
  integer vector of scanID will be generated as \code{1:nSamples}.  This
  ID is used to index scans in the netCDF file.

  The map file will be converted to a
  \code{\link{SnpAnnotationDataFrame}}.  This SNP annotation will
  include the definitions of A and B alleles in the netCDF file (either
  as provided or determined from the data as described above).  A unique
  integer snpID will be generated for each SNP, which is used to index
  SNPs in the netCDF file.
  
  Note that the default values of \code{ncdfXYchromCode=24}, \code{ncdfYchromCode=25},
  and \code{ncdfUchromCode=27} correspond to the default chromosome
  codes for \code{\link{NcdfGenotypeReader}} and
  \code{\link{SnpAnnotationDataFrame}}, and  are different from the values used by
  PLINK (Y=24, XY=25, U=0).  If the netCDF file is created with
  different chromosome codes by specifying these arguments, one must
  also specify the chromosome codes when opening the file,
  e.g. \code{NcdfGenotypeReader(ncdfFile, XYchromCode=25,
  YchromCode=24)}.
  
  \code{nSamples} is used to allocate space in the netCDF file.  A
  warning will be issued if the number of lines read in the ped file is
  different from this number.
}
\references{
Please see
  \url{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped} for
  more information on PLINK files. 
}
\author{Stephanie Gogarten}
\seealso{\code{\link{plinkWrite}}, \code{\link{plinkCheck}}
}
\examples{
library(GWASdata)
pedfile <- system.file("extdata", "illumina_subj.ped", package="GWASdata")
mapfile <- system.file("extdata", "illumina_subj.map", package="GWASdata")
ncfile <- tempfile()
scanfile <- tempfile()
snpfile <- tempfile()
plinkToNcdf(pedfile, mapfile, nSamples=43, ncdfFile=ncfile,
   snpAnnotFile=snpfile, scanAnnotFile=scanfile)

nc <- NcdfGenotypeReader(ncfile)
scanAnnot <- getobj(scanfile)
snpAnnot <- getobj(snpfile)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
prefix <- sub(".ped", "", pedfile, fixed=TRUE)
log <- tempfile()
stopifnot(plinkCheck(genoData, prefix, log, alleleA.col="alleleA",
  alleleB.col="alleleB"))
close(genoData)

# provide allele coding with extended map file
# .bim might have SNPs in different order than .map
bimfile <- system.file("extdata", "illumina_subj.bim", package="GWASdata")
bim <- read.table(bimfile, as.is=TRUE, header=FALSE)
map <- read.table(mapfile, as.is=TRUE, header=FALSE)
snp.match <- match(map[,2], bim[,2])
map <- cbind(map, bim[snp.match, 5:6])
mapfile.ext <- tempfile()
write.table(map, file=mapfile.ext, quote=FALSE, row.names=FALSE, col.names=FALSE)
# use chromosome codes that match PLINK
plinkToNcdf(pedfile, mapfile, nSamples=43, ncdfFile=ncfile,
   snpAnnotFile=snpfile, scanAnnotFile=scanfile,
   ncdfYchromCode=24, ncdfXYchromCode=25)

# must specify different chromosome codes in NcdfGenotypeReader
# appending "L" ensures the codes are integers, as required
nc <- NcdfGenotypeReader(ncfile, YchromCode=24L, XYchromCode=25L)
scanAnnot <- getobj(scanfile)
snpAnnot <- getobj(snpfile)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
stopifnot(plinkCheck(genoData, prefix, log, alleleA.col="alleleA",
  alleleB.col="alleleB"))
close(genoData)

file.remove(ncfile, scanfile, snpfile, log, mapfile.ext)
}
\keyword{manip}
