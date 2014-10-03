\name{IntensityData-class}
\docType{class}

\alias{IntensityData-class}
\alias{IntensityData}
\alias{show,IntensityData-method}
\alias{open,IntensityData-method}
\alias{close,IntensityData-method}
\alias{nsnp,IntensityData-method}
\alias{nscan,IntensityData-method}
\alias{hasSnpAnnotation,IntensityData-method}
\alias{hasScanAnnotation,IntensityData-method}
\alias{getSnpID,IntensityData-method}
\alias{getChromosome,IntensityData-method}
\alias{getPosition,IntensityData-method}
\alias{getScanID,IntensityData-method}
\alias{hasSex,IntensityData-method}
\alias{getSex,IntensityData-method}
\alias{hasSnpVariable,IntensityData-method}
\alias{getSnpVariable,IntensityData-method}
\alias{getSnpVariableNames,IntensityData-method}
\alias{hasScanVariable,IntensityData-method}
\alias{getScanVariable,IntensityData-method}
\alias{getScanVariableNames,IntensityData-method}
\alias{hasVariable,IntensityData-method}
\alias{getVariable,IntensityData-method}
\alias{getQuality,IntensityData-method}
\alias{getX,IntensityData-method}
\alias{getY,IntensityData-method}
\alias{getBAlleleFreq,IntensityData-method}
\alias{getLogRRatio,IntensityData-method}
\alias{hasQuality,IntensityData-method}
\alias{hasX,IntensityData-method}
\alias{hasY,IntensityData-method}
\alias{hasBAlleleFreq,IntensityData-method}
\alias{hasLogRRatio,IntensityData-method}
\alias{autosomeCode,IntensityData-method}
\alias{XchromCode,IntensityData-method}
\alias{XYchromCode,IntensityData-method}
\alias{YchromCode,IntensityData-method}
\alias{MchromCode,IntensityData-method}

\title{Class IntensityData}

\description{
  The IntensityData class is a container for storing intensity data from a
  genome-wide association study together with the metadata associated
  with the subjects and SNPs involved in the study.
}

\details{
  The IntensityData class consists of three slots: data, snp annotation,
  and scan annotation.  There may be
  multiple scans associated with a subject (e.g. duplicate scans for
  quality control), hence the use of "scan" as one dimension of the
  data.
  Snp and scan annotation are optional, but if included in the
  IntensityData object, their unique integer ids (snpID and scanID) are
  checked against the ids stored in the data file to ensure consistency.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{IntensityData(data, snpAnnot=NULL, scanAnnot=NULL)}:

      \code{data} must be a \code{\link{GdsIntensityReader}} or \code{\link{NcdfIntensityReader}} object.

      \code{snpAnnot}, if not \code{NULL}, must be a
      \code{\link{SnpAnnotationDataFrame}} or \code{\link{SnpAnnotationSQLite}} object.

      \code{scanAnnot}, if not \code{NULL}, must be a
      \code{\link{ScanAnnotationDataFrame}} or \code{\link{ScanAnnotationSQLite}} object.

      The \code{IntensityData} constructor creates and returns a
      IntensityData instance, ensuring that data, snpAnnot, and scanAnnot
      are internally consistent.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is an IntensityData object.
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
      \code{getQuality(object, snp, scan)}: Extracts quality scores.  
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{getX(object, snp, scan)}: Extracts X intensity values.  
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{getY(object, snp, scan)}: Extracts Y intensity values.  
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{getBAlleleFreq(object, snp, scan)}: Extracts B allele
      frequency values.  
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{getLogRRatio(object, snp, scan)}: Extracts Log R Ratio values.  
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
  \code{\link{ScanAnnotationDataFrame}},
  \code{\link{GdsIntensityReader}},
  \code{\link{NcdfIntensityReader}},
  \code{\link{GenotypeData}}
}

\examples{
library(GWASdata)
file <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
gds <- GdsIntensityReader(file)

# object without annotation
intenData <- IntensityData(gds)

# object with annotation
data(illuminaSnpADF, illuminaScanADF)
intenData <-  IntensityData(gds, snpAnnot=illuminaSnpADF, scanAnnot=illuminaScanADF)

# dimensions
nsnp(intenData)
nscan(intenData)

# get snpID and chromosome
snpID <- getSnpID(intenData)
chrom <- getChromosome(intenData)

# get positions only for chromosome 22
pos22 <- getPosition(intenData, index=(chrom == 22))

# get other annotations
if (hasSex(intenData)) sex <- getSex(intenData)
plate <- getScanVariable(intenData, "plate")
rsID <- getSnpVariable(intenData, "rsID")

# get all snps for first scan
x <- getX(intenData, snp=c(1,-1), scan=c(1,1))

# starting at snp 100, get 10 snps for the first 5 scans
x <- getX(intenData, snp=c(100,10), scan=c(1,5))

close(intenData)
}

\keyword{methods}
\keyword{classes}

