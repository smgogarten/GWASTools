\name{gdsImputedDosage}
\alias{gdsImputedDosage}
\alias{gdsCheckImputedDosage}
\title{Create and check a GDS file with imputed dosages}
\description{
  These functions create or check a GDS file and corresponding
  annotation for imputed dosages from IMPUTE2, BEAGLE, or MaCH.
}
\usage{
gdsImputedDosage(input.files, gds.filename, chromosome,
                  input.type=c("IMPUTE2", "BEAGLE", "MaCH"), 
                  input.dosage=FALSE, block.size=5000,
                  snp.annot.filename="dosage.snp.RData",
                  scan.annot.filename="dosage.scan.RData",
                  verbose=TRUE, zipflag="ZIP.max", 
                  genotypeDim="snp,scan",
                  scan.df=NULL, snp.exclude=NULL, snp.id.start=1)
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
  \item{gds.filename}{Character string with name of output GDS
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
  \item{genotypeDim}{character string specifying genotype dimensions of gds file. Either "snp,scan" or "scan,snp"}
  \item{zipflag}{compression method for GDS nodes, can be "", "ZIP", "ZIP.fast", "ZIP.default", "ZIP.max"}
  \item{scan.df}{data frame specifying which samples to include in the output GDS files, with optional scanIDs already assigned. See details.}
  \item{snp.exclude}{vector of integers specifying which SNPs to exclude from the GDS file.}
  \item{snp.id.start}{Starting index for snpID.}
  \item{genoData}{A \code{\link{GenotypeData}} object from a GDS file created with \code{\link{gdsImputedDosage}}.}
  \item{snpAnnot}{The \code{\link{SnpAnnotationDataFrame}} created by \code{\link{gdsImputedDosage}}}
  \item{scanAnnot}{The \code{\link{ScanAnnotationDataFrame}} created by \code{\link{gdsImputedDosage}}}  
  \item{tolerance}{tolerance for checking differences against input files}

}
\details{
  Input files can contain either imputed dosages or genotype
  probabilities, specified by the \code{input.dosage} flag.
  In either case, the GDS file will store dosage of
  the A allele in the "genotype" variable.  All SNPs are assumed to be
  on the same chromosome, which is indicated by the \code{chromosome}
  argument.

  The \code{scan.df} argument allows the user to specify what samples should be included in the GDS files and an optional sampleID-scanID mapping. \code{scan.df} is a data frame with required column \code{sampleID}. The function attempts to match the given sampleID in the \code{scan.df} data frame with a unique sampleID in the input files. The format of sampleID is different for different input types:
  \itemize{
  \item{IMPUTE2:}{ "ID_1 ID_2" as given in the sample file, where IDs are separated by a space}
  \item{BEAGLE:}{ Column header names corresponding to that sample in .dose or .gprobs file}
  \item{MaCH:}{ The first column of the .mlprob or .mlprob file}
  }
  
  The \code{snp.names} argument allows the user to specify the which SNPs should be included in the GDS files. However, \code{snp.names} must be in the same order as SNPs occur in the imputation files; this option therefore only allows selection of SNPs, not reordering of SNPs. The ordering is checked and an error is thrown if the SNP names are not in order, but due to the design of imputation files, this may not occur until well into the GDS file population. The user can specify the starting snpID by setting \code{snp.id.start}, and included SNPs are numbered sequentially starting with \code{snp.id.start}. For IMPUTE2 data, snp.names must correspond to the second column of the .gprobs file.

  Minimal SNP and scan annotation are created from the input files and stored in
  RData format in \code{snp.annot.filename} and \code{scan.annot.filename}.

  Currently supported input file types are IMPUTE2, BEAGLE, and MaCH.
}
\references{
  IMPUTE2: \url{http://mathgen.stats.ox.ac.uk/impute/impute_v2.html}

  BEAGLE:
  \url{http://faculty.washington.edu/browning/beagle/beagle.html}

  MaCH: \url{http://www.sph.umich.edu/csg/abecasis/MACH/tour/imputation.html}
}
\author{Adrienne Stilp, Stephanie Gogarten}
\seealso{\code{GdsGenotypeReader}, \code{\link{GenotypeData}}, \code{\link{assocTestRegression}}}
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

geno <- getGenotype(genoData)
getVariable(genoData, "alleleA")
getVariable(genoData, "alleleB")

# association test with imputed dosages
scanAnnot$status <- sample(0:1, nrow(scanAnnot), replace=TRUE)
genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
assoc <- assocTestRegression(genoData, outcome="status", model.type="logistic",
                             gene.action.list="additive", dosage=TRUE)
head(assoc)
close(genoData)


# BEAGLE - genotype probabilities
probfile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.gprobs",
                      package="GWASdata")
markfile <- system.file("extdata", "imputation", "BEAGLE", "hapmap.markers",
                    package="GWASdata")
gdsImputedDosage(input.files=c(probfile, markfile), gds.filename=gdsfile, chromosome=22,
                  input.type="BEAGLE", input.dosage=FALSE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)

# BEAGLE - dosage
dosefile <- system.file("extdata", "imputation", "BEAGLE", "example.hapmap.unphased.bgl.dose",
                    package="GWASdata")
gdsImputedDosage(input.files=c(dosefile, markfile), gds.filename=gdsfile, chromosome=22,
                  input.type="BEAGLE", input.dosage=TRUE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)


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

# MaCH - dosage
dosefile <- system.file("extdata", "imputation", "MaCH", "mach1.out.mldose",
                        package="GWASdata")
gdsImputedDosage(input.files=c(dosefile, markfile, posfile), gds.filename=gdsfile, chromosome=22,
                  input.type="MaCH", input.dosage=TRUE,
                  snp.annot.filename=snpfile, scan.annot.filename=scanfile)

unlink(c(gdsfile, snpfile, scanfile))
}
\keyword{manip}