\name{ncdfAddData}
\alias{ncdfAddData}
\alias{ncdfAddIntensity}
\alias{ncdfCheckGenotype}
\alias{ncdfCheckIntensity}

\title{
Write genotypic calls and/or associated metrics to a netCDF file
}
\description{
Genotypic calls and/or associated quantitative variables (e.g. quality score, intensities) are read from text files and written to an existing netCDF file in which those variables were defined previously.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usage{
ncdfAddData(path = ".", ncdf.filename, 
            snp.annotation, scan.annotation,
            sep.type, skip.num, col.total, col.nums, 
            scan.name.in.file, scan.start.index = 1,
            diagnostics.filename = "ncdfAddData.diagnostics.RData",
            verbose = TRUE)

ncdfAddIntensity(path = ".",  ncdf.filename,
                 snp.annotation, scan.annotation, 
                 scan.start.index = 1, n.consecutive.scans = -1,  
                 diagnostics.filename = "ncdfAddIntensity.diagnostics.RData",
                 verbose = TRUE)

ncdfCheckGenotype(path = ".", ncdf.filename, 
                  snp.annotation, scan.annotation, 
                  sep.type, skip.num, col.total, col.nums, 
                  scan.name.in.file, check.scan.index, n.scans.loaded,
                  diagnostics.filename = "ncdfCheckGenotype.diagnostics.RData",
                  verbose = TRUE)

ncdfCheckIntensity(path = ".", intenpath = ".", ncdf.filename, 
                   snp.annotation, scan.annotation, 
                   sep.type, skip.num, col.total, col.nums, 
                   scan.name.in.file, check.scan.index,
                   n.scans.loaded, affy.inten = FALSE,
                   diagnostics.filename = "ncdfCheckIntensity.diagnostics.RData",
                   verbose = TRUE)
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{path}{Path to the raw text files.
}

  \item{intenpath}{Path to the raw text files containing intensity, if
    "inten.file" is given in scan.annotation.
}

  \item{ncdf.filename}{Name of the netCDF file in which to write the data.
}

  \item{snp.annotation}{SNP annotation data.frame containing
    SNPs in the same order as those in the snp dimension of the netCDF
    file.  Column names must be "snpID" (integer ID) and "snpName",
    where snpName matches the snp ids inside the raw genoypic data files.
}

\item{scan.annotation}{Scan annotation data.frame with columns "scanID"
  (integer id of scan in the netCDF file), "scanName", (sample name
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
    c("snp", "sample", "geno", "a1", "a2", "qs", "x", "y", "rawx",
    "rawy", "r", "theta", "ballelefreq", "logrratio"). The element "snp"
    is the column of SNP ids, "sample" is sample ids, "geno" is diploid
    genotype (in AB format), "a1" and "a2" are alleles 1 and 2 (in AB
    format), "qs" is quality score, "x" and "y" are normalized
    intensities, "rawx" and "rawy" are raw intensities, "r" is the sum
    of normalized intensities, "theta" is angular polar coordinate,
    "ballelefreq" is the B allele frequency, and "logrratio" is the Log
    R Ratio.
}

  \item{scan.name.in.file}{An indicator for the presence of sample name within the file. A value of 1 indicates a column with repeated values of the sample name (Illumina format), -1 indicates sample name embedded in a column heading (Affymetrix format) and 0 indicates no sample name inside the raw data file.
}

  \item{scan.start.index}{A numeric value containing the index of the sample dimension of the netCDF file at which to begin writing.
  }

  \item{n.consecutive.scans}{The number of consecutive "sampleID" indices for which to write intensity values, beginning at scan.start.index (which equals the number of "ALLELE_SUMMARY" files to process). When n.consecutive.scans=-1, all samples from scan.start.index to the total number will be processed.}
  
  \item{check.scan.index}{An integer vector containing the indices of
  the sample dimension of the netCDF file to check.
}

  \item{n.scans.loaded}{Number of scans loaded in the netCDF file.}

  \item{affy.inten}{Logical value indicating whether Affy intensities
  are in separate files from quality scores.  If \code{TRUE}, must also specify intenpath.}
  
  \item{diagnostics.filename}{Name of the output file to save diagnostics.}

  \item{verbose}{Logical value specifying whether to show progress information.}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{
These functions read genotypic and associated data from raw text
files. The files to be read and processed are specified in the sample
annotation. \code{ncdfAddData} expects one file per sample, with each
file having one row of data per SNP probe. The \code{col.nums} argument
allows the user to select and identify specific fields for writing to
the netCDF file. Illumina text files and Affymetrix ".CHP" files can be
used here (but not Affymetrix "ALLELE_SUMMARY" files).

A SNP annotation data.frame is a pre-requisite for this function. It has the same number of rows (one per SNP) as the raw text file and a column of SNP names matching those within the raw text file. It also has a column of integer SNP ids matching the values (in order) of the "snp" dimension of the netCDF file.

A sample annotation data.frame is also a pre-requisite. It has one row per sample with columns corresponding to sample name (as it occurs within the raw text file), name of the raw text file for that sample and an integer sample id (to be written as the "sampleID" variable in the netCDF file).

The genotype calls in the raw text file may be either one column of diploid calls or two columns of allele calls. The function takes calls in AB format and converts them to a numeric code indicating the number of "A" alleles in the genotype (i.e. AA=2, AB=1, BB=0 and missing=-1).

While each raw text file is being read, the functions check for errors and irregularities and records the results in a list of vectors. If any problem is detected, that raw text file is skipped.

\code{ncdfAddIntensity} uses \code{scan.start.index} and \code{n.consecutive.scans} to identify the set of integer sample ids for input (from the netCDF file). It then uses the sample annotation data.frame to identify the corresponding sample names and "ALLELE_SUMMARY" file names to read.
The "ALLELE_SUMMARY" files have two rows per SNP, one for X (A allele) and one for Y (B allele). These are reformatted to one row per SNP and and ordered according to the SNP integer id in the netCDF file. The correspondence between SNP names in the "ALLELE_SUMMARY" file and the SNP integer ids is made using the SNP annotation data.frame.

\code{ncdfCheckGenotype} and \code{ncdfCheckIntensity} check the
contents of netCDF files against raw text files.

These functions use the \pkg{\link{ncdf}} library, which provides an interface between R and netCDF.


}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
The netCDF file specified in argument \code{ncdf.filename} is populated with genotype calls and/or associated quantitative variables. 
A list of diagnostics with the following components is returned. Each vector has one element per raw text file processed.
\item{read.file}{ A vector indicating whether (1) or not (0) each file was read successfully.}
\item{row.num}{	A vector of the number of rows read from each file. These should all be the same and equal to the number of rows in the SNP annotation data.frame.}
\item{samples}{ A list of vectors containing the unique sample names in the sample column of each raw text file. Each vector should have just one element.}
\item{sample.match}{ A vector indicating whether (1) or not (0) the sample name inside the raw text file matches that in the sample annotation data.frame}
\item{missg}{ A list of vectors containing the unique character string(s) for missing genotypes (i.e. not AA,AB or BB) for each raw text file.}
\item{snp.chk}{ A vector indicating whether (1) or not (0) the raw text file has the expected set of SNP names (i.e. matching those in the SNP annotation data.frame).}
\item{chk}{	A vector indicating whether (1) or not (0) all previous
  checks were successful and the data were written to the netCDF file.}

\code{ncdfCheckGenotypes} returns the following additional list items.
\item{snp.order}{A vector indicating whether (1) or not (0) the snp ids
  are in the same order in each file.}
\item{geno.chk}{A vector indicating whether (1) or not (0) the genotypes
  in the netCDF match the text file.}

\code{ncdfCheckIntensity} returns the following additional list items.
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\author{Cathy Laurie
}

\note{These functions were modeled after similar code written by Thomas Lumley.
}


\seealso{\pkg{\link{ncdf}}, \code{\link{ncdfCreate}}, \code{\link{ncdfSubset}}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
library(GWASdata)

#############
# Illumina - genotype file
#############
# first create empty netCDF
data(illumina_snp_annot)
snpAnnot <- illumina_snp_annot
data(illumina_scan_annot)
scanAnnot <- illumina_scan_annot[1:3,] # subset of samples for testing
ncfile <- tempfile()
ncdfCreate(snpAnnot, ncfile, variables="genotype",
                n.samples=nrow(scanAnnot))

# add data
path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
snpAnnot <- snpAnnot[,c("snpID", "rsID")]
names(snpAnnot) <- c("snpID", "snpName")
scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "file")]
names(scanAnnot) <- c("scanID", "scanName", "file")
col.nums <- as.integer(c(1,2,12,13))
names(col.nums) <- c("snp", "sample", "a1", "a2")
diagfile <- tempfile()
res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type=",",
                     skip.num=11, col.total=21, col.nums=col.nums,
                     scan.name.in.file=1, diagnostics.filename=diagfile)

file.remove(diagfile)
file.remove(ncfile)

#############
# Affymetrix - genotype file
#############
# first create empty netCDF
data(affy_snp_annot)
snpAnnot <- affy_snp_annot
data(affy_scan_annot)
scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
ncfile <- tempfile()
ncdfCreate(snpAnnot, ncfile, variables="genotype",
                n.samples=nrow(scanAnnot))

# add data
path <- system.file("extdata", "affy_raw_data", package="GWASdata")
snpAnnot <- snpAnnot[,c("snpID", "probeID")]
names(snpAnnot) <- c("snpID", "snpName")
scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "chpFile")]
names(scanAnnot) <- c("scanID", "scanName", "file")
col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
diagfile <- tempfile()
res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                     skip.num=1, col.total=6, col.nums=col.nums,
                     scan.name.in.file=-1, diagnostics.filename=diagfile)
file.remove(diagfile)

# check
diagfile <- tempfile()
res <- ncdfCheckGenotype(path, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                       skip.num=1, col.total=6, col.nums=col.nums,
                       scan.name.in.file=-1, check.scan.index=1:3,
                       n.scans.loaded=3, diagnostics.filename=diagfile)
file.remove(diagfile)
file.remove(ncfile)

#############
# Affymetrix - intensity file
#############
# first create empty netCDF
snpAnnot <- affy_snp_annot
scanAnnot <- affy_scan_annot[1:3,] # subset of samples for testing
ncfile <- tempfile()
ncdfCreate(snpAnnot, ncfile, variables=c("quality","X","Y"),
                n.samples=nrow(scanAnnot))

# add sampleID and quality
path <- system.file("extdata", "affy_raw_data", package="GWASdata")
snpAnnot <- snpAnnot[,c("snpID", "probeID")]
names(snpAnnot) <- c("snpID", "snpName")
scanAnnot1 <- scanAnnot[,c("scanID", "genoRunID", "chpFile")]
names(scanAnnot1) <- c("scanID", "scanName", "file")
col.nums <- as.integer(c(2,4)); names(col.nums) <- c("snp", "qs")
diagfile <- tempfile()
res <- ncdfAddData(path, ncfile, snpAnnot, scanAnnot1, sep.type="\t",
                     skip.num=1, col.total=6, col.nums=col.nums,
                     scan.name.in.file=-1, diagnostics.filename=diagfile)
file.remove(diagfile)

# add intensity
scanAnnot2 <- scanAnnot[,c("scanID", "genoRunID", "alleleFile")]
names(scanAnnot2) <- c("scanID", "scanName", "file")
diagfile <- tempfile()
res <- ncdfAddIntensity(path, ncfile, snpAnnot, scanAnnot2,
                        diagnostics.filename=diagfile)
file.remove(diagfile)

# check
intenpath <- system.file("extdata", "affy_raw_data", package="GWASdata")
scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "chpFile", "alleleFile")]
names(scanAnnot) <- c("scanID", "scanName", "file", "inten.file")
diagfile <- tempfile()
res <- ncdfCheckIntensity(path, intenpath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                        skip.num=1, col.total=6, col.nums=col.nums,
                        scan.name.in.file=-1, check.scan.index=1:3,
                        n.scans.loaded=3, affy.inten=TRUE,
                        diagnostics.filename=diagfile)
  
file.remove(diagfile)
file.remove(ncfile)
}
\keyword{manip}
