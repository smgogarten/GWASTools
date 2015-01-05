\name{NcdfIntensityReader}
\docType{class}

\alias{NcdfIntensityReader-class}
\alias{NcdfIntensityReader}
\alias{getVariable,NcdfIntensityReader-method}
\alias{getSnpID,NcdfIntensityReader-method}
\alias{getChromosome,NcdfIntensityReader-method}
\alias{getPosition,NcdfIntensityReader-method}
\alias{getScanID,NcdfIntensityReader-method}
\alias{getQuality,NcdfIntensityReader-method}
\alias{getX,NcdfIntensityReader-method}
\alias{getY,NcdfIntensityReader-method}
\alias{getBAlleleFreq,NcdfIntensityReader-method}
\alias{getLogRRatio,NcdfIntensityReader-method}
\alias{hasQuality,NcdfIntensityReader-method}
\alias{hasX,NcdfIntensityReader-method}
\alias{hasY,NcdfIntensityReader-method}
\alias{hasBAlleleFreq,NcdfIntensityReader-method}
\alias{hasLogRRatio,NcdfIntensityReader-method}
\alias{nsnp,NcdfIntensityReader-method}
\alias{nscan,NcdfIntensityReader-method}
\alias{autosomeCode,NcdfIntensityReader-method}
\alias{XchromCode,NcdfIntensityReader-method}
\alias{XYchromCode,NcdfIntensityReader-method}
\alias{YchromCode,NcdfIntensityReader-method}
\alias{MchromCode,NcdfIntensityReader-method}

\title{Class NcdfIntensityReader}

\description{
  The NcdfIntensityReader class is an extension of the NcdfReader class
  specific to reading intensity data stored in NetCDF files.
}

\section{Extends}{
  \code{\link{NcdfReader}}
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{NcdfIntensityReader(filename)}:

      \code{filename} must be the path to a NetCDF file.  The NetCDF
      file must contain the following variables:
      \itemize{
	\item 'snp': a coordinate variable with a unique
	integer vector of snp ids
	\item 'chromosome':  integer
	chromosome values of dimension 'snp'
	\item 'position': integer position values of dimension 'snp'
	\item 'sampleID': a unique integer vector of scan ids with
	dimension 'sample'
      }
      Default values for chromosome codes are 1-22=autosome, 23=X, 24=XY, 25=Y,
      26=M.  The defaults may be changed with the arguments \code{autosomeCode},
      \code{XchromCode}, \code{XYchromCode}, \code{YchromCode}, and
      \code{MchromCode}.

      The NetCDF file should also contain at least one of the following
      variables with dimensions ('snp','sample'):
      \itemize{
	\item 'quality': quality score
	\item 'X': X intensity
	\item 'Y': Y intensity
	\item 'BAlleleFreq': B allele frequency
	\item 'LogRRatio': Log R Ratio
      }

      The \code{NcdfIntensityReader} constructor creates and returns a
      NcdfIntensityReader instance pointing to this file.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a NcdfIntensityReader object.
  \code{snp} and
      \code{scan} indicate which elements to return along the snp and
      scan dimensions.  They must be integer vectors of the form (start,
      count), where start is the index of the first data element to read
      and count is the number of elements to read.  A value of '-1' for
      count indicates that the entire dimension should be read.  If snp
      and/or is scan omitted, the entire variable is read.

      See \code{\link{NcdfReader}} for additional methods.

  \describe{
    \item{}{
      \code{nsnp(object)}: The number of SNPs in the NetCDF file.
    }
    \item{}{
      \code{nscan(object)}: The number of scans in the NetCDF file.
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
      (1:22,X,XY,Y,M,U).  "U" stands for "Unknown" and is the value
      given to any chromosome code not falling in the other categories.
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
      \code{getQuality(object, snp, scan)}: Extracts quality scores.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{hasQuality(object)}:
      Returns \code{TRUE} if the GDS file contains a variable 'quality'.
    }
    \item{}{
      \code{getX(object, snp, scan)}: Extracts X intensity.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{hasX(object)}:
      Returns \code{TRUE} if the GDS file contains a variable 'X'.
    }
    \item{}{
      \code{getY(object, snp, scan)}: Extracts Y intensity.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{hasY(object)}:
      Returns \code{TRUE} if the GDS file contains a variable 'Y'.
    }
    \item{}{
      \code{getBAlleleFreq(object, snp, scan)}: Extracts B allele frequency.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{hasBAlleleFreq(object)}:
      Returns \code{TRUE} if the GDS file contains a variable 'BAlleleFreq'.
    }
    \item{}{
      \code{getLogRRatio(object, snp, scan)}: Extracts Log R Ratio.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.
    }
    \item{}{
      \code{hasLogRRatio(object)}:
      Returns \code{TRUE} if the GDS file contains a variable 'LogRRatio'.
    }
    \item{}{
      \code{getVariable(object, varname, snp, scan)}: Returns the
      contents of the variable \code{varname}.
      The result is a vector or matrix, depending on the number
      of dimensions in the returned values.  Missing values are
      represented as \code{NA}.  If the variable is not found in the NetCDF
      file, returns \code{NULL}.
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

\seealso{\code{\link{NcdfReader}},
  \code{\link{NcdfGenotypeReader}},
  \code{\link{GenotypeData}}, \code{\link{IntensityData}}
}

\examples{
file <- system.file("extdata", "illumina_qxy.nc", package="GWASdata")
nc <- NcdfIntensityReader(file)

# dimensions
nsnp(nc)
nscan(nc)

# get snpID and chromosome
snpID <- getSnpID(nc)
chrom <- getChromosome(nc)

# get positions only for chromosome 22
pos22 <- getPosition(nc, index=(chrom == 22))

# get all snps for first scan
x <- getX(nc, snp=c(1,-1), scan=c(1,1))

# starting at snp 100, get 10 snps for the first 5 scans
x <- getX(nc, snp=c(100,10), scan=c(1,5))

close(nc)
}

\keyword{methods}
\keyword{classes}

