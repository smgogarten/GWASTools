\name{GWASTools-package}
\alias{GWASTools-package}
\alias{GWASTools}
\docType{package}
\title{Tools for Genome Wide Association Studies}

\description{
  This package contains tools for facilitating cleaning (quality control
  and quality assurance) and analysis of GWAS data.
}

\details{
  GWASTools provides a set of classes for storing data and annotation
  from Genome Wide Association studies, and a set of functions for data
  cleaning and analysis that operate on those classes.

  Genotype and intensity data are stored in NetCDF files, so it is
  possible to analyze data sets that are too large to be contained in
  memory.  The \code{\link{NcdfReader}} class provides a generic
  interface to the NetCDF files (utilizing the \code{\link{ncdf}} package), and the
  \code{\link{NcdfGenotypeReader}} and \code{\link{NcdfIntensityReader}}
  classes provide specific methods to access genotype and intensity
  data.

  Two sets of classes for annotation are provided.  
  \code{\link{SnpAnnotationDataFrame}} and 
  \code{\link{ScanAnnotationDataFrame}} extend
  \code{\link{AnnotatedDataFrame}} and provide in-memory containers for SNP and
  scan annotation and metadata.    
  \code{\link{SnpAnnotationSQLite}} and 
  \code{\link{ScanAnnotationSQLite}} provide interfaces to SNP and scan
  annotation and metadata stored in SQLite databases.
  
  The \code{\link{GenotypeData}} and \code{\link{IntensityData}} classes
  combine genotype or intensity data with SNP and scan annotation,
  ensuring that the data in the NetCDF files is consistent with
  annotation through unique SNP and scan IDs.  A majority of the
  functions in the GWASTools package take \code{\link{GenotypeData}}
  and/or \code{\link{IntensityData}} objects as arguments.
}

\author{
  Stephanie Gogarten, Cathy Laurie, Tushar Bhangale, Matt Conomos, Cecilia Laurie, Caitlin McHugh, Ian Painter, Xiuwen Zheng, Rohit Swarnkar
  
  Maintainer: Stephanie Gogarten \email{sdmorris@u.washington.edu}
}

\references{
Laurie, C. C., Doheny, K. F., Mirel, D. B., Pugh, E. W., Bierut, L. J.,
Bhangale, T., Boehm, F., Caporaso, N. E., Cornelis, M. C., Edenberg,
H. J., Gabriel, S. B., Harris, E. L., Hu, F. B., Jacobs, K. B., Kraft,
P., Landi, M. T., Lumley, T., Manolio, T. A., McHugh, C., Painter, I.,
Paschall, J., Rice, J. P., Rice, K. M., Zheng, X., and Weir, B. S., for
the GENEVA Investigators (2010), Quality control and quality assurance
in genotypic data for genome-wide association studies. Genetic
Epidemiology, 34: 591-602. doi: 10.1002/gepi.20516
}

\keyword{package}
