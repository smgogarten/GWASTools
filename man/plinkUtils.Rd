\name{plinkUtils}
\alias{plinkWrite}
\alias{plinkCheck}
\title{Utilities to create and check PLINK files}
\description{
\code{plinkWrite} creates ped and map format files (used by PLINK) from a \code{\link{GenotypeData}}
object.
\code{plinkCheck} checks whether a set of ped and map files has
identical data to a \code{\link{GenotypeData}} object.
}
\usage{
plinkWrite(genoData, pedFile="testPlink", family.col="family",
  individual.col="scanID", father.col="father", mother.col="mother",
  phenotype.col=NULL, alleleA.col=NULL, alleleB.col=NULL,
  rs.col="rsID", mapdist.col=NULL, scan.exclude=NULL,
  blockSize=100)

plinkCheck(genoData, pedFile, family.col="family",
  individual.col="scanID", father.col="father", mother.col="mother",
  phenotype.col=NULL, alleleA.col=NULL, alleleB.col=NULL,
  rs.col="rsID")
}
\arguments{
  \item{genoData}{A \code{\link{GenotypeData}} object with scan and SNP annotation.}
  \item{pedFile}{prefix for PLINK files (pedFile.ped, pedFile.map)}
  \item{family.col}{name of the column in the scan annotation that contains family ID of the sample}
  \item{individual.col}{name of the column in the scan annotation that contains individual ID of the sample}
  \item{father.col}{name of the column in the scan annotation that contains father ID of the sample}
  \item{mother.col}{name of the column in the scan annotation that contains mother ID of the sample}
  \item{phenotype.col}{name of the column in the scan annotation that
    contains phenotype variable (e.g. case control statue) of the sample}
  \item{alleleA.col}{name of the column in the SNP annototation that
  contains allele A for the SNP}
  \item{alleleB.col}{name of the column in the SNP annototation that
  contains allele B for the SNP}
  \item{rs.col}{name of the column in the SNP annotation that contains rs ID (or some other ID) for the SNP}
  \item{mapdist.col}{name of the column in the SNP annotation that
    contains genetic distance in Morgans for the SNP}
  \item{scan.exclude}{vector of scanIDs to exclude from PLINK file}
  \item{blockSize}{Number of samples to read from \code{genoData} at a time}
}
\details{
If \code{alleleA.col} or \code{alleleB.col} is NULL, genotypes are
written as "A A", "A B", "B B" (or "0 0" for missing data).

If \code{phenotype.col=NULL}, \code{plinkWrite} will use "-9"
for writing phenotype data and \code{plinkCheck} will omit checking this
column.

If \code{mapdist.col=NULL}, \code{plinkWrite} will use "0"
for writing this column in the map file and \code{plinkCheck} will omit checking this
column.

\code{plinkCheck} first reads the map file and checks for SNP
mismatches (chromosome, rsID, and/or position).  If any mismatches are
found, a file "snp_mismatch.log" is written.
If no SNP errors were found, the code reads the ped file line by line and stops at the first
mismatch.  SNPs and sample order is not required to be the same as in \code{genoData}.

These utilities convert between chromosome coding in
\code{\link{GenotypeData}}, which by default is 24=XY, 25=Y, and PLINK
chromosome coding, which is 24=Y, 25=X.

Larger blockSize will improve speed but will require more RAM.
}
\value{
  \code{plinkCheck} returns \code{TRUE} if the PLINK files contain
  identical data to \code{genoData}, and \code{FALSE} is a mismatch is encountered.
}
\references{
Please see
  \url{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped} for
  more information on the ped and map files. 
}
\author{Stephanie Gogarten, Tushar Bhangale}
\examples{
library(GWASdata)
ncfile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
data(illumina_snp_annot)
data(illumina_scan_annot)

# family data is required
illumina_scan_annot$family <- 0
illumina_scan_annot$father <- 0
illumina_scan_annot$mother <- 0

genoData <- GenotypeData(NcdfGenotypeReader(ncfile),
  scanAnnot=ScanAnnotationDataFrame(illumina_scan_annot),
  snpAnnot=SnpAnnotationDataFrame(illumina_snp_annot))

pedfile <- tempfile()
plinkWrite(genoData, pedfile,
  alleleA.col="alleleA", alleleB.col="alleleB")
plinkCheck(genoData, pedfile,
  alleleA.col="alleleA", alleleB.col="alleleB")

# use default AB allele coding and exclude samples
plinkWrite(genoData, pedfile, scan.exclude=c(281, 283),
  blockSize=10)
plinkCheck(genoData, pedfile)
#samples not found in Ped:
#281
#283

unlink(paste(pedfile, "*", sep="."))
}
\keyword{manip}
