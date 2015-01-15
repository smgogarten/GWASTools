\name{GenotypeData-class}
\docType{class}

\alias{GenotypeData-class}
\alias{GenotypeData}
\alias{show,GenotypeData-method}
\alias{open,GenotypeData-method}
\alias{close,GenotypeData-method}
\alias{nsnp,GenotypeData-method}
\alias{nscan,GenotypeData-method}
\alias{hasSnpAnnotation,GenotypeData-method}
\alias{hasScanAnnotation,GenotypeData-method}
\alias{getSnpID,GenotypeData-method}
\alias{getChromosome,GenotypeData-method}
\alias{getPosition,GenotypeData-method}
\alias{getAlleleA,GenotypeData-method}
\alias{getAlleleB,GenotypeData-method}
\alias{getScanID,GenotypeData-method}
\alias{hasSex,GenotypeData-method}
\alias{getSex,GenotypeData-method}
\alias{hasSnpVariable,GenotypeData-method}
\alias{getSnpVariable,GenotypeData-method}
\alias{getSnpVariableNames,GenotypeData-method}
\alias{hasScanVariable,GenotypeData-method}
\alias{getScanVariable,GenotypeData-method}
\alias{getScanVariableNames,GenotypeData-method}
\alias{hasVariable,GenotypeData-method}
\alias{getVariable,GenotypeData-method}
\alias{getGenotype,GenotypeData-method}
\alias{getGenotypeSelection,GenotypeData-method}
\alias{autosomeCode,GenotypeData-method}
\alias{XchromCode,GenotypeData-method}
\alias{XYchromCode,GenotypeData-method}
\alias{YchromCode,GenotypeData-method}
\alias{MchromCode,GenotypeData-method}

\title{Class GenotypeData}

\description{
  The GenotypeData class is a container for storing genotype data from a
  genome-wide association study together with the metadata associated
  with the subjects and SNPs involved in the study.
}

\details{
  The GenotypeData class consists of three slots: data, snp annotation,
  and scan annotation.  There may be
  multiple scans associated with a subject (e.g. duplicate scans for
  quality control), hence the use of "scan" as one dimension of the
  data.
  Snp and scan annotation are optional, but if included in the
  GenotypeData object, their unique integer ids (snpID and scanID) are
  checked against the ids stored in the data slot to ensure consistency.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenotypeData(data, snpAnnot=NULL, scanAnnot=NULL)}:

      \code{data} must be an \code{\link{NcdfGenotypeReader}}, \code{\link{GdsGenotypeReader}}, or \code{\link{MatrixGenotypeReader}} object.

      \code{snpAnnot}, if not \code{NULL}, must be a
      \code{\link{SnpAnnotationDataFrame}} or \code{\link{SnpAnnotationSQLite}} object.

      \code{scanAnnot}, if not \code{NULL}, must be a
      \code{\link{ScanAnnotationDataFrame}} or \code{\link{ScanAnnotationSQLite}} object.

      The \code{GenotypeData} constructor creates and returns a
      GenotypeData instance, ensuring that data, snpAnnot, and scanAnnot
      are internally consistent.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a GenotypeData object.

  \describe{
    \item{}{
      \code{nsnp(object)}: The number of SNPs in the data.
    }
    \item{}{
      \code{nscan(object)}: The number of scans in the data.
    }
    \item{}{
      \code{getSnpID(object, index)}: A unique integer vector of snp
      IDs.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.
    }
    \item{}{
      \code{getChromosome(object, index, char=FALSE)}: A vector of
      chromosomes.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.
      If \code{char=FALSE} (default), returns an integer vector.
      If \code{char=TRUE}, returns a character vector with elements in
      (1:22,X,XY,Y,M,U).
    }
    \item{}{
      \code{getPosition(object, index)}: An integer vector of base pair
      positions.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.
    }
    \item{}{
      \code{getAlleleA(object, index)}: A character vector of A alleles.
      The optional \code{index} is a logical or
      integer vector specifying elements to extract.
    }
    \item{}{
      \code{getAlleleB(object, index)}: A character vector of B alleles.
      The optional \code{index} is a logical or
      integer vector specifying elements to extract.
    }
    \item{}{
      \code{getScanID(object, index)}: A unique integer vector of scan
      IDs.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.
    }
    \item{}{
      \code{getSex(object, index)}: A character vector of sex, with values 'M'
      or 'F'.    The optional \code{index} is a logical or
      integer vector specifying elements to extract.
    }
    \item{}{
      \code{hasSex(object)}: Returns \code{TRUE} if the column 'sex' is present in
      \code{object}.
    }
    \item{}{
      \code{getGenotype(object, snp=c(1,-1), scan=c(1,-1), char=FALSE, sort=TRUE, drop=TRUE, ...)}: Extracts genotype values
      (number of A alleles).
      \code{snp} and
      \code{scan} indicate which elements to return along the snp and
      scan dimensions.  They must be integer vectors of the form (start,
      count), where start is the index of the first data element to read
      and count is the number of elements to read.  A value of '-1' for
      count indicates that the entire dimension should be read.
      If \code{drop=TRUE} the result is coerced to the lowest possible dimension.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values and the value of \code{drop}.  
      Missing values are
      represented as \code{NA}.  If \code{char=TRUE}, genotypes are
      returned as characters of the form "A/B".  If \code{sort=TRUE},
      alleles are lexographically sorted ("G/T" instead of "T/G").
    }
    \item{}{
      \code{getGenotypeSelection(object, snp=NULL, scan=NULL, char=FALSE, sort=TRUE, drop=TRUE, ...)}: May
    be used only if the data slot contains a
    \code{\link{GdsGenotypeReader}} or \code{\link{MatrixGenotypeReader}} object.
      Extracts genotype values
      (number of A alleles).
      \code{snp} and
      \code{scan} may be integer or logical vectors indicating which elements to return along the snp and
      scan dimensions.  Unlike \code{getGenotype}, the values requested need not be in contiguous blocks.
      Other arguments are identical to \code{getGenotype}.
    }
    \item{}{
      \code{getSnpVariable(object, varname, index)}: Returns the snp
      annotation variable \code{varname}.
      The optional \code{index} is a logical or
      integer vector specifying elements to extract.
    }
    \item{}{
      \code{getSnpVariableNames(object)}: Returns a character vector
      with the names of the columns in the snp annotation.
    }
    \item{}{
      \code{hasSnpVariable(object, varname)}: Returns \code{TRUE} if the
      variable \code{varname} is present in the snp annotation.
    }
    \item{}{
      \code{getScanVariable(object, varname, index)}: Returns the scan
      annotation variable \code{varname}.
      The optional \code{index} is a logical or
      integer vector specifying elements to extract.
    }
    \item{}{
      \code{getScanVariableNames(object)}: Returns a character vector
      with the names of the columns in the scan annotation.
    }
    \item{}{
      \code{hasScanVariable(object, varname)}: Returns \code{TRUE} if the
      variable \code{varname} is present in the scan annotation.
    }
    \item{}{
      \code{getVariable(object, varname, drop=TRUE, ...)}: Extracts the
      contents of the variable \code{varname} from the data.
      If \code{drop=TRUE} the result is coerced to the lowest possible dimension.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values and the value of \code{drop}.  Missing values are
      represented as \code{NA}.  If the variable is not found, returns \code{NULL}.
    }
    \item{}{
      \code{hasVariable(object, varname)}: Returns \code{TRUE} if
      the data contains contains \code{varname}, \code{FALSE} if not.
    }
    \item{}{
      \code{hasSnpAnnotation(object)}: Returns \code{TRUE} if the snp
      annotation slot is not \code{NULL}.
    }
    \item{}{
      \code{hasScanAnnotation(object)}: Returns \code{TRUE} if the scan
      annotation slot is not \code{NULL}.
    }
    \item{}{
      \code{open(object)}: Opens a connection to the data.
    }
    \item{}{
      \code{close(object)}: Closes the data connection.
    }
     \item{}{
       \code{autosomeCode(object)}: Returns the integer codes for the
       autosomes.
     }
     \item{}{
       \code{XchromCode(object)}: Returns the integer code for the X
       chromosome.
     }
     \item{}{
       \code{XYchromCode(object)}: Returns the integer code for the
       pseudoautosomal region.
     }
     \item{}{
       \code{YchromCode(object)}: Returns the integer code for the Y
       chromosome.
     }
     \item{}{
       \code{MchromCode(object)}: Returns the integer code for
       mitochondrial SNPs.
     }
  }
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SnpAnnotationDataFrame}},
  \code{\link{SnpAnnotationSQLite}},
  \code{\link{ScanAnnotationDataFrame}},
  \code{\link{ScanAnnotationSQLite}},
  \code{\link{GdsGenotypeReader}},
  \code{\link{NcdfGenotypeReader}},
  \code{\link{MatrixGenotypeReader}},
  \code{\link{IntensityData}}
}

\examples{
library(GWASdata)
file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(file)

# object without annotation
genoData <- GenotypeData(gds)

# object with annotation
data(illuminaSnpADF)
data(illuminaScanADF)
# need to rebuild old SNP annotation object to get allele methods
snpAnnot <- SnpAnnotationDataFrame(pData(illuminaSnpADF))
genoData <-  GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot=illuminaScanADF)

# dimensions
nsnp(genoData)
nscan(genoData)

# get snpID and chromosome
snpID <- getSnpID(genoData)
chrom <- getChromosome(genoData)

# get positions only for chromosome 22
pos22 <- getPosition(genoData, index=(chrom == 22))

# get other annotations
if (hasSex(genoData)) sex <- getSex(genoData)
plate <- getScanVariable(genoData, "plate")
rsID <- getSnpVariable(genoData, "rsID")

# get all snps for first scan
geno <- getGenotype(genoData, snp=c(1,-1), scan=c(1,1))

# starting at snp 100, get 10 snps for the first 5 scans
geno <- getGenotype(genoData, snp=c(100,10), scan=c(1,5))
geno

# return genotypes as "A/B" rather than number of A alleles
geno <- getGenotype(genoData, snp=c(100,10), scan=c(1,5), char=TRUE)
geno

close(genoData)

#--------------------------------------
# An example using a non-human organism
#--------------------------------------
# Chicken has 38 autosomes, Z, and W.  Male is ZZ, female is ZW.
# Define sex chromosomes as X=Z and Y=W.
ncfile <- tempfile()
simulateGenotypeMatrix(n.snps=10, n.chromosomes=40, n.samples=5,
                       ncdf.filename=ncfile)
nc <- NcdfGenotypeReader(ncfile, autosomeCode=1:38L,
                         XchromCode=39L, YchromCode=40L,
                         XYchromCode=41L, MchromCode=42L)
table(getChromosome(nc))
table(getChromosome(nc, char=TRUE))

# SNP annotation
snpdf <- data.frame(snpID=getSnpID(nc),
                    chromosome=getChromosome(nc),
                    position=getPosition(nc))
snpAnnot <- SnpAnnotationDataFrame(snpdf, autosomeCode=1:38L,
                         XchromCode=39L, YchromCode=40L,
                         XYchromCode=41L, MchromCode=42L)
varMetadata(snpAnnot)[,"labelDescription"] <-
  c("unique integer ID",
    "chromosome coded as 1:38=autosomes, 39=Z, 40=W",
    "base position")

# reverse sex coding to get proper counting of sex chromosome SNPs
scandf <- data.frame(scanID=1:5, sex=c("M","M","F","F","F"),
                     stringsAsFactors=FALSE)
scanAnnot <- ScanAnnotationDataFrame(scandf)
varMetadata(scanAnnot)[,"labelDescription"] <-
  c("unique integer ID",
    "sex coded as M=female and F=male")

genoData <- GenotypeData(nc, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
afreq <- alleleFrequency(genoData)
# frequency of Z chromosome in females ("M") and males ("F")
afreq[snpAnnot$chromosome == 39, c("M","F")]
# frequency of W chromosome in females ("M") and males ("F")
afreq[snpAnnot$chromosome == 40, c("M","F")]

close(genoData)
unlink(ncfile)
}

\keyword{methods}
\keyword{classes}

