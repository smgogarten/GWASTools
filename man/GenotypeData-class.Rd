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

      \code{data} must be an \code{\link{NcdfGenotypeReader}} or \code{\link{MatrixGenotypeReader}} object.

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
  \code{snp} and
      \code{scan} indicate which elements to return along the snp and
      scan dimensions.  They must be integer vectors of the form (start,
      count), where start is the index of the first data element to read
      and count is the number of elements to read.  A value of '-1' for
      count indicates that the entire dimension should be read.  If snp
      and/or is scan omitted, the entire variable is read.
      
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
      \code{getGenotype(object, snp, scan)}: Extracts genotype values
      (number of A alleles).  
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
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
      \code{getVariable(object, varname, snp, scan)}: Extracts the
      contents of the variable \code{varname} from the data.  
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
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
  \code{\link{NcdfReader}},
  \code{\link{NcdfGenotypeReader}},
  \code{\link{MatrixGenotypeReader}},
  \code{\link{IntensityData}}
}

\examples{
library(GWASdata)
file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
nc <- NcdfGenotypeReader(file)

# object without annotation
genoData <- GenotypeData(nc)

# object with annotation
data(affySnpADF)
data(affyScanADF)
genoData <-  GenotypeData(nc, snpAnnot=affySnpADF, scanAnnot=affyScanADF)

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

close(genoData)
}

\keyword{methods}
\keyword{classes}

