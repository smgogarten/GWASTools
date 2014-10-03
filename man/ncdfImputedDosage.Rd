\name{ncdfImputedDosage}
\alias{ncdfImputedDosage}
\title{Create a NetCDF file with imputed dosages}
\description{
  This function creates a NetCDF file and corresponding
  annotation for imputed dosages from IMPUTE2, BEAGLE, or MaCH.
}
\usage{
ncdfImputedDosage(input.files, ncdf.filename, chromosome,
                  input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                  input.dosage=FALSE, block.size=5000,
                  snp.annot.filename="dosage.snp.RData",
                  scan.annot.filename="dosage.scan.RData",
                  verbose=TRUE)
}
\arguments{
  \item{input.files}{A character vector of input files.  The first file
    should always be genotypes (either probabilities or dosages).  Files
    for each input type should be as follows:
    \itemize{
      \item IMPUTE2: 1) .gens, 2) .samples
      \item BEAGLE: 1) .grobs or .dose, 2) .markers
      \item MaCH: 1) .mlprob or .mldose, 2) .mlinfo, 3) file with
      columns named "SNP" and "position" giving base pair position of
      all SNPs
    }
  }
  \item{ncdf.filename}{Character string with name of output NetCDF
    file.}
  \item{chromosome}{Chromosome corresponding to the SNPs in the genotype
    file.  Character codes will be mapped to integer values as follows:
    "X"->23, "XY"->24, "Y"-> 25, "M","MT"->26.}
  \item{input.type}{Format of input files.  Accepted file types are
    "IMPUTE2", "BEAGLE", and "MaCH".}
  \item{input.dosage}{Logical for whether the genotype file
    (\code{input.files[1]}) contains dosages.  If \code{FALSE} (default), the
    genotype file is assumed to contain genotype probabilities.}
  \item{block.size}{Number of lines to read at once.}
  \item{snp.annot.filename}{Output .RData file for storing a \code{\link{SnpAnnotationDataFrame}}.}
  \item{scan.annot.filename}{Output .RData file for storing a
    \code{\link{ScanAnnotationDataFrame}}.}
  \item{verbose}{Logical for whether to print progress messages.}
}
\details{
  Input files can contain either imputed dosages or genotype
  probabilities, specified by the \code{input.dosage} flag.
  In either case, the NetCDF file will store dosage of
  the A allele in the "genotype" variable.  All SNPs are assumed to be
  on the same chromosome, which is indicated by the \code{chromosome}
  argument.

  SNP and scan annotation are created from the input files and stored in
  RData format in \code{snp.annot.filename} and \code{scan.annot.filename}.

  Currently supported input file types are IMPUTE2, BEAGLE, and MaCH.
}
\references{
  IMPUTE2: \url{http://mathgen.stats.ox.ac.uk/impute/impute_v2.html}

  BEAGLE:
  \url{http://faculty.washington.edu/browning/beagle/beagle.html}

  MaCH: \url{http://www.sph.umich.edu/csg/abecasis/MACH/tour/imputation.html}
}
\author{Stephanie Gogarten}
\seealso{\code{\link{createDataFile}}, \code{\link{NcdfGenotypeReader}}, \code{\link{GenotypeData}}, \code{\link{assocTestRegression}}}
\examples{
ncfile <- tempfile()
snpfile <- tempfile()
scanfile <- tempfile()

# IMPUTE2
probfile <- system.file("extdata", "imputation", "IMPUTE2", "example.chr22.study.gens",
                        package="GWASdata")
sampfile <- system.file("extdata", "imputation", "IMPUTE2", "example.study.samples",
                        package="GWASdata")
ncdfImputedDosage(input.files=c(probfile, sampfile), ncdf.filename=ncfile, chromosome=22,
                  input.type="IMPUTE2", input.dosage=FALSE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)
nc <- NcdfGenotypeReader(ncfile)
scanAnnot <- getobj(scanfile)
snpAnnot <- getobj(snpfile)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
geno <- getGenotype(genoData)
getVariable(genoData, "alleleA")
getVariable(genoData, "alleleB")

# association test with imputed dosages
scanAnnot$status <- sample(0:1, nrow(scanAnnot), replace=TRUE)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
assoc <- assocTestRegression(genoData, outcome="status", model.type="logistic",
                             gene.action.list="additive", dosage=TRUE)
head(assoc)
close(genoData)


# BEAGLE - genotype probabilities
probfile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.gprobs",
                      package="GWASdata")
markfile <- system.file("extdata", "imputation", "BEAGLE", "hapmap.markers",
                    package="GWASdata")
ncdfImputedDosage(input.files=c(probfile, markfile), ncdf.filename=ncfile, chromosome=22,
                  input.type="BEAGLE", input.dosage=FALSE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)

# BEAGLE - dosage
dosefile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.dose",
                    package="GWASdata")
ncdfImputedDosage(input.files=c(dosefile, markfile), ncdf.filename=ncfile, chromosome=22,
                  input.type="BEAGLE", input.dosage=TRUE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)


# MaCH - genotype probabilities
probfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlprob",
                        package="GWASdata")
markfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlinfo",
                        package="GWASdata")
posfile <- system.file("extdata", "imputation", "MaCH", "mach1.snp.position",
                        package="GWASdata")
ncdfImputedDosage(input.files=c(probfile, markfile, posfile), ncdf.filename=ncfile, chromosome=22,
                  input.type="MaCH", input.dosage=FALSE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)

# MaCH - dosage
dosefile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mldose",
                        package="GWASdata")
ncdfImputedDosage(input.files=c(dosefile, markfile, posfile), ncdf.filename=ncfile, chromosome=22,
                  input.type="MaCH", input.dosage=TRUE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)

unlink(c(ncfile, snpfile, scanfile))
}
\keyword{manip}