\name{createDataFile}
\alias{createDataFile}
%\alias{ncdfAddIntensity}
\alias{checkGenotypeFile}
\alias{checkIntensityFile}


\title{
Write genotypic calls and/or associated metrics to a GDS or netCDF file.
}

\description{
Genotypic calls and/or associated quantitative variables (e.g. quality score, intensities) are read from text files and written to a GDS or netCDF file.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\usage{
createDataFile(path = ".", filename, file.type=c("gds", "ncdf"),
               variables="genotype", snp.annotation, scan.annotation,
               sep.type, skip.num, col.total, col.nums,
               scan.name.in.file, scan.start.index = 1,
	       precision="single", compress="ZIP.max",
	       array.name = NULL, genome.build = NULL,
               diagnostics.filename = "createDataFile.diagnostics.RData",
               verbose = TRUE)

%ncdfAddIntensity(path = ".",  filename, file.type=c("gds", "ncdf"),
%                 snp.annotation, scan.annotation,
%                 scan.start.index = 1, n.consecutive.scans = -1,
%                 diagnostics.filename = "ncdfAddIntensity.diagnostics.RData",
%                 verbose = TRUE)

checkGenotypeFile(path = ".", filename, file.type=c("gds", "ncdf"),
                  snp.annotation, scan.annotation,
                  sep.type, skip.num, col.total, col.nums,
                  scan.name.in.file, check.scan.index, n.scans.loaded,
                  diagnostics.filename = "checkGenotypeFile.diagnostics.RData",
                  verbose = TRUE)

checkIntensityFile(path = ".", intenpath = ".", filename, file.type=c("gds", "ncdf"),
                   snp.annotation, scan.annotation,
                   sep.type, skip.num, col.total, col.nums,
                   scan.name.in.file, check.scan.index,
                   n.scans.loaded, affy.inten = FALSE,
                   diagnostics.filename = "checkIntensityFile.diagnostics.RData",
                   verbose = TRUE)
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\arguments{

  \item{path}{Path to the raw text files.
}

  \item{intenpath}{Path to the raw text files containing intensity, if
    "inten.file" is given in scan.annotation.
}

\item{filename}{The name of the genotype GDS or netCDF file to create}

\item{file.type}{The type of file to create ("gds" or "ncdf")}

\item{variables}{A character vector containing the names of the
  variables to create (must be one or more of \code{c("genotype", "quality", "X", "Y", "rawX", "rawY", "R", "Theta", "BAlleleFreq", "LogRRatio")})}

\item{snp.annotation}{Snp annotation dataframe with columns "snpID",
  "chromosome",  "position" and "snpName".  snpID should be a unique integer
  vector, sorted with
  respect to chromosome and position.
  snpName should match the snp identifiers inside the raw genoypic data files
  If \code{file.type="gds"},
  optional columns  "alleleA", and "alleleB" will be
  written if present.}

\item{scan.annotation}{Scan annotation data.frame with columns "scanID"
  (integer id of genotyping instance), "scanName", (sample name
  inside the raw data file) and "file" (corresponding raw data file name).
}

  \item{sep.type}{Field separator in the raw text files.
}

  \item{skip.num}{Number of rows to skip, which should be all rows preceding the genotypic or quantitative data (including the header).
}

  \item{col.total}{Total number of columns in the raw text files.
}

\item{col.nums}{An integer vector indicating which columns of the raw text
    file contain variables for input.  \code{names(col.nums)} must be a subset of
    c("snp", "sample", "geno", "a1", "a2", "quality", "X", "Y", "rawX",
    "rawY", "R", "Theta", "BAlleleFreq", "LogRRatio"). The element "snp"
    is the column of SNP ids, "sample" is sample ids, "geno" is diploid
    genotype (in AB format), "a1" and "a2" are alleles 1 and 2 (in AB
    format), "quality" is quality score, "X" and "Y" are normalized
    intensities, "rawX" and "rawY" are raw intensities, "R" is the sum
    of normalized intensities, "Theta" is angular polar coordinate,
    "BAlleleFreq" is the B allele frequency, and "LogRRatio" is the Log
    R Ratio.
}

  \item{scan.name.in.file}{An indicator for the presence of sample name within the file. A value of 1 indicates a column with repeated values of the sample name (Illumina format), -1 indicates sample name embedded in a column heading (Affymetrix format) and 0 indicates no sample name inside the raw data file.
}

  \item{scan.start.index}{A numeric value containing the index of the sample dimension of the GDS or netCDF file at which to begin writing.
  }

%  \item{n.consecutive.scans}{The number of consecutive "sampleID" indices for which to write intensity values, beginning at scan.start.index (which equals the number of "ALLELE_SUMMARY" files to process). When n.consecutive.scans=-1, all samples from scan.start.index to the total number will be processed.}

  \item{check.scan.index}{An integer vector containing the indices of
  the sample dimension of the GDS or netCDF file to check.}

  \item{n.scans.loaded}{Number of scans loaded in the GDS or netCDF file.}

  \item{affy.inten}{Logical value indicating whether Affy intensities
  are in separate files from quality scores.  If \code{TRUE}, must also specify intenpath.}

\item{precision}{A character value indicating whether floating point
  numbers should be stored as "double" or "single" precision.}

\item{compress}{The compression level for variables in a GDS file (see
\code{\link{add.gdsn}} for options.}

\item{array.name}{Name of the array, to be stored as an attribute in the
  netCDF file.}

\item{genome.build}{Genome build used in determining chromosome and
  position, to be stored as an attribute in the netCDF file.}

  \item{diagnostics.filename}{Name of the output file to save diagnostics.}

  \item{verbose}{Logical value specifying whether to show progress information.}


}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\details{
These functions read genotypic and associated data from raw text
files. The files to be read and processed are specified in the sample
annotation. \code{createDataFile} expects one file per sample, with each
file having one row of data per SNP probe. The \code{col.nums} argument
allows the user to select and identify specific fields for writing to
the GDS or netCDF file. Illumina text files and Affymetrix ".CHP" files can be
used here (but not Affymetrix "ALLELE_SUMMARY" files).

A SNP annotation data.frame is a pre-requisite for this function. It has the same number of rows (one per SNP) as the raw text file and a column of SNP names matching those within the raw text file. It also has a column of integer SNP ids to be used as a unique key for each SNP in the GDS or netCDF file.

A sample annotation data.frame is also a pre-requisite. It has one row per sample with columns corresponding to sample name (as it occurs within the raw text file), name of the raw text file for that sample and an integer sample id (to be written as the "sampleID" variable in the GDS or netCDF file).

The genotype calls in the raw text file may be either one column of diploid calls or two columns of allele calls. The function takes calls in AB format and converts them to a numeric code indicating the number of "A" alleles in the genotype (i.e. AA=2, AB=1, BB=0 and missing=-1).

While each raw text file is being read, the functions check for errors and irregularities and records the results in a list of vectors. If any problem is detected, that raw text file is skipped.

%\code{ncdfAddIntensity} uses \code{scan.start.index} and \code{n.consecutive.scans} to identify the set of integer sample ids for input (from the netCDF file). It then uses the sample annotation data.frame to identify the corresponding sample names and "ALLELE_SUMMARY" file names to read.
%The "ALLELE_SUMMARY" files have two rows per SNP, one for X (A allele) and one for Y (B allele). These are reformatted to one row per SNP and and ordered according to the SNP integer id in the netCDF file. The correspondence between SNP names in the "ALLELE_SUMMARY" file and the SNP integer ids is made using the SNP annotation data.frame.

\code{checkGenotypeFile} and \code{checkIntensityFile} check the
contents of GDS or netCDF files against raw text files.

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
The GDS or netCDF file specified in argument \code{filename} is populated with genotype calls and/or associated quantitative variables.
A list of diagnostics with the following components is returned. Each vector has one element per raw text file processed.
\item{read.file}{ A vector indicating whether (1) or not (0) each file was read successfully.}
\item{row.num}{	A vector of the number of rows read from each file. These should all be the same and equal to the number of rows in the SNP annotation data.frame.}
\item{samples}{ A list of vectors containing the unique sample names in the sample column of each raw text file. Each vector should have just one element.}
\item{sample.match}{ A vector indicating whether (1) or not (0) the sample name inside the raw text file matches that in the sample annotation data.frame}
\item{missg}{ A list of vectors containing the unique character string(s) for missing genotypes (i.e. not AA,AB or BB) for each raw text file.}
\item{snp.chk}{ A vector indicating whether (1) or not (0) the raw text file has the expected set of SNP names (i.e. matching those in the SNP annotation data.frame).}
\item{chk}{	A vector indicating whether (1) or not (0) all previous
  checks were successful and the data were written to the netCDF file.}

\code{checkGenotypeFile} returns the following additional list items.
\item{snp.order}{A vector indicating whether (1) or not (0) the snp ids
  are in the same order in each file.}
\item{geno.chk}{A vector indicating whether (1) or not (0) the genotypes
  in the netCDF match the text file.}

\code{checkIntensityFile} returns the following additional list items.
\item{qs.chk}{A vector indicating whether (1) or not (0) the quality scores
  in the netCDF match the text file.}
\item{read.file.inten}{ A vector indicating whether (1) or not (0) each
  intensity file was read successfully (if intensity files are separate).}
\item{sample.match.inten}{ A vector indicating whether (1) or not (0)
  the sample name inside the raw text file matches that in the sample
  annotation data.frame (if intensity files are separate). }
\item{rows.equal}{A vector indicating whether (1) or not (0) the number
  of rows read from each file are the same and equal to the number of rows in
  the SNP annotation data.frame (if intensity files are separate).}
\item{snp.chk.inten}{ A vector indicating whether (1) or not (0) the raw text
  file has the expected set of SNP names (i.e. matching those in the SNP
  annotation data.frame) (if intensity files are separate).}
\item{inten.chk}{A vector for each intensity variable indicating whether (1) or not (0) the
  intensities in the netCDF match the text file.}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\author{
Stephanie Gogarten, Cathy Laurie
}
\seealso{\pkg{\link{gdsfmt}}, \pkg{\link{ncdf}}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
library(GWASdata)

#############
# Illumina - genotype file
#############
ncfile <- tempfile()
path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
data(illumina_snp_annot, illumina_scan_annot)
snpAnnot <- illumina_snp_annot[,c("snpID", "rsID", "chromosome",
	                          "position", "alleleA", "alleleB")]
names(snpAnnot)[2] <-  "snpName"
# subset of samples for testing
scanAnnot <- illumina_scan_annot[1:3, c("scanID", "genoRunID", "file")]
names(scanAnnot)[2] <- "scanName"
col.nums <- as.integer(c(1,2,12,13))
names(col.nums) <- c("snp", "sample", "a1", "a2")
diagfile <- tempfile()
res <- createDataFile(path, ncfile, file.type="ncdf", variables="genotype",
                     snpAnnot, scanAnnot, sep.type=",",
                     skip.num=11, col.total=21, col.nums=col.nums,
                     scan.name.in.file=1, diagnostics.filename=diagfile)

file.remove(diagfile)
file.remove(ncfile)

#############
# Affymetrix - genotype file
#############
ncfile <- tempfile()
path <- system.file("extdata", "affy_raw_data", package="GWASdata")
data(affy_snp_annot, affy_scan_annot)
snpAnnot <- affy_snp_annot[,c("snpID", "probeID", "chromosome", "position")]
names(snpAnnot)[2] <- "snpName"
# subset of samples for testing
scanAnnot <- affy_scan_annot[1:3, c("scanID", "genoRunID", "chpFile")]
names(scanAnnot)[2:3] <- c("scanName", "file")
col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
diagfile <- tempfile()
res <- createDataFile(path, ncfile, file.type="ncdf", variables="genotype",
                     snpAnnot, scanAnnot, sep.type="\t",
                     skip.num=1, col.total=6, col.nums=col.nums,
                     scan.name.in.file=-1, diagnostics.filename=diagfile)
file.remove(diagfile)

# check
diagfile <- tempfile()
res <- checkGenotypeFile(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                       skip.num=1, col.total=6, col.nums=col.nums,
                       scan.name.in.file=-1, check.scan.index=1:3,
                       n.scans.loaded=3, diagnostics.filename=diagfile)
file.remove(diagfile)
file.remove(ncfile)
}

\keyword{manip}
