\name{gdsCheckImputedDosage}
\alias{gdsCheckImputedDosage}
\title{Checks an imputed dosage GDS file against source files}
\description{
  This function checks a GDS file created with gdsImputedDosage against the source files.
  IMPUTE2, BEAGLE, and MaCH file formats are supported.s
}
\usage{
gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot, 
                      input.files, chromosome,
                      input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                      input.dosage=FALSE, block.size=5000,
                      verbose=TRUE, 
                      snp.exclude=NULL,
                      snp.id.start=1,
                      tolerance=1e-4)
}
\arguments{
  \item{genoData}{A \code{\link{GenotypeData}} object from a GDS file created with \code{\link{gdsImputedDosage}}.}
  \item{snpAnnot}{The \code{\link{SnpAnnotationDataFrame}} created by \code{\link{gdsImputedDosage}}}
  \item{scanAnnot}{The \code{\link{ScanAnnotationDataFrame}} created by \code{\link{gdsImputedDosage}}}  
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
  \item{chromosome}{Chromosome corresponding to the SNPs in the genotype
    file.  Character codes will be mapped to integer values as follows:
    "X"->23, "XY"->24, "Y"-> 25, "M","MT"->26.}
  \item{input.type}{Format of input files.  Accepted file types are
    "IMPUTE2", "BEAGLE", and "MaCH".}
  \item{input.dosage}{Logical for whether the genotype file
    (\code{input.files[1]}) contains dosages.  If \code{FALSE} (default), the
    genotype file is assumed to contain genotype probabilities.}
  \item{block.size}{Number of lines to read at once.}
  \item{verbose}{Logical for whether to print progress messages.}
  \item{snp.exclude}{vector of integers specifying which SNPs to exclude from the GDS file.}
  \item{snp.id.start}{Starting index for snpID.}
  \item{tolerance}{tolerance for checking differences against input files}
}
\details{
  See \code{\link{gdsImputedDosage}} for more detailed description of \code{input.files}, \code{input.type}, \code{input.dosage}, \code{snp.exclude}, and \code{snp.id.start}.
}
\references{
  IMPUTE2: \url{http://mathgen.stats.ox.ac.uk/impute/impute_v2.html}

  BEAGLE:
  \url{http://faculty.washington.edu/browning/beagle/beagle.html}

  MaCH: \url{http://www.sph.umich.edu/csg/abecasis/MACH/tour/imputation.html}
}
\author{Adrienne Stilp}
\seealso{\code{gdsImputedDosage}, \code{GdsGenotypeReader}, \code{\link{GenotypeData}}}
\examples{
gdsfile <- tempfile()
snpfile <- tempfile()
scanfile <- tempfile()

# IMPUTE2
probfile <- system.file("extdata", "imputation", "IMPUTE2", "example.chr22.study.gens",
                        package="GWASdata")
sampfile <- system.file("extdata", "imputation", "IMPUTE2", "example.study.samples",
                        package="GWASdata")
gdsImputedDosage(input.files=c(probfile, sampfile), gds.filename=gdsfile, chromosome=22,
                  input.type="IMPUTE2", input.dosage=FALSE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)
gds <- GdsGenotypeReader(gdsfile)
scanAnnot <- getobj(scanfile)
snpAnnot <- getobj(snpfile)
genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot,
                      input.files=c(probfile, sampfile), chromosome=22,
                      input.type="IMPUTE2", input.dosage=FALSE)



# BEAGLE - genotype probabilities
probfile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.gprobs",
                      package="GWASdata")
markfile <- system.file("extdata", "imputation", "BEAGLE", "hapmap.markers",
                    package="GWASdata")
gdsImputedDosage(input.files=c(probfile, markfile), gds.filename=gdsfile, chromosome=22,
                  input.type="BEAGLE", input.dosage=FALSE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)
gds <- GdsGenotypeReader(gdsfile)
scanAnnot <- getobj(scanfile)
snpAnnot <- getobj(snpfile)
genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot,
                      input.files=c(probfile, markfile), chromosome=22,
                      input.type="BEAGLE", input.dosage=FALSE)
close(genoData)

# BEAGLE - dosage
dosefile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.dose",
                    package="GWASdata")
gdsImputedDosage(input.files=c(dosefile, markfile), gds.filename=gdsfile, chromosome=22,
                  input.type="BEAGLE", input.dosage=TRUE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)
gds <- GdsGenotypeReader(gdsfile)
scanAnnot <- getobj(scanfile)
snpAnnot <- getobj(snpfile)
genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot,
                      input.files=c(dosefile, markfile), chromosome=22,
                      input.type="BEAGLE", input.dosage=TRUE)
close(genoData)


# MaCH - genotype probabilities
probfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlprob",
                        package="GWASdata")
markfile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mlinfo",
                        package="GWASdata")
posfile <- system.file("extdata", "imputation", "MaCH", "mach1.snp.position",
                        package="GWASdata")
gdsImputedDosage(input.files=c(probfile, markfile, posfile), gds.filename=gdsfile, chromosome=22,
                  input.type="MaCH", input.dosage=FALSE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)
gds <- GdsGenotypeReader(gdsfile)
scanAnnot <- getobj(scanfile)
snpAnnot <- getobj(snpfile)
genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot,
                      input.files=c(probfile, markfile, posfile), chromosome=22,
                      input.type="MaCH", input.dosage=FALSE)
close(genoData)


# MaCH - dosage
dosefile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mldose",
                        package="GWASdata")
gdsImputedDosage(input.files=c(dosefile, markfile, posfile), gds.filename=gdsfile, chromosome=22,
                  input.type="MaCH", input.dosage=TRUE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)
gds <- GdsGenotypeReader(gdsfile)
scanAnnot <- getobj(scanfile)
snpAnnot <- getobj(snpfile)
genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

gdsCheckImputedDosage(genoData, snpAnnot, scanAnnot,
                      input.files=c(dosefile, markfile, posfile), chromosome=22,
                      input.type="MaCH", input.dosage=TRUE)
close(genoData)

unlink(c(gdsfile, snpfile, scanfile))
}
\keyword{manip}