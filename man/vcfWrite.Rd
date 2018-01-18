\name{vcfWrite}
\alias{vcfWrite}
\alias{vcfCheck}
\alias{genoDataAsVCF}
\title{Utility to write VCF file}
\description{
\code{genoDataAsVCF} creates a \code{\link[VariantAnnotation]{VCF-class}} object.

\code{vcfWrite} writes a VCF file from a \code{\link{GenotypeData}}
object.

\code{vcfCheck} compares the genotypes in a VCF file to
the corresponding genotypes in \code{genoData}.

}
\usage{
genoDataAsVCF(genoData, sample.col="scanID",
              id.col="snpID", qual.col=NULL, filter.cols=NULL,
              info.cols=NULL, scan.exclude=NULL, snp.exclude=NULL,
              scan.order=NULL, ref.allele=NULL)

vcfWrite(genoData, vcf.file="out.vcf", sample.col="scanID",
         id.col="snpID", qual.col=NULL, filter.cols=NULL,
         info.cols=NULL, scan.exclude=NULL, snp.exclude=NULL,
         scan.order=NULL, ref.allele=NULL, block.size=1000, verbose=TRUE)

vcfCheck(genoData, vcf.file,
         sample.col="scanID", id.col="snpID",
         scan.exclude=NULL, snp.exclude=NULL,
         block.size=1000, verbose=TRUE)
}
\arguments{
  \item{genoData}{A \code{\link{GenotypeData}} object with scan and SNP annotation.}
  \item{vcf.file}{Filename for the output VCF file.}
  \item{sample.col}{name of the column in the scan annotation to use as sample IDs in the VCF file}
  \item{id.col}{name of the column in the SNP annotation to use as "ID" column in the VCF file}
  \item{qual.col}{name of the column in the SNP annotation to use as "QUAL" column in the VCF file}
  \item{filter.cols}{vector of column names in the SNP annotation to use as "FILTER" column in the VCF file.  These columns should be logical vectors, with \code{TRUE} for SNPs to be filtered.  Any SNPs with a value of \code{FALSE} for all filter columns will be set to "PASS".}
  \item{info.cols}{vector of column names in the SNP annotation to concatenate for the "INFO" column in the VCF file.}
  \item{scan.exclude}{vector of scanIDs to exclude from creation or
         checking of VCF file}
  \item{snp.exclude}{vector of snpIDs to exclude from creation or
         checking of VCF file}
  \item{scan.order}{vector of scanIDs to include in VCF file, in the order in which they should be written}
  \item{ref.allele}{vector of "A" or "B" values indicating where allele A or allele B should be the reference allele for each SNP.  Default is to use allele A as the reference allele.}
  \item{block.size}{Number of SNPs to read from \code{genoData} at a
    time}
  \item{verbose}{logical for whether to show progress information.}
}
\details{
	REF will be alleleA and ALT will be alleleB.

	vcfCheck compares the genotypes (diploid only) in a VCF file to
	the corresponding genotypes in \code{genoData}.  It stops with
	an error when it detects a discordant genotype.  It assumes that
	the "ID" column of the VCF file has unique values that can be
	matched with a column in the SNP annotation. The VCF file may
	contain additional samples or SNPs not present in the
	\code{genoData}; these records will be automaticlaly excluded from the check.
	Users may exclude additional SNPs and samples (i.e. those
	overlapping with \code{genoData}) using the \code{scan.exclude}
	and \code{snp.exclude} arguments.
}
\references{
	The variant call format and VCFtools.
	Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE,
	Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group.
	Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.
}
\author{Stephanie Gogarten, Michael Lawrence, Sarah Nelson}
\seealso{\code{\link[SNPRelate]{snpgdsVCF2GDS}}}
\examples{
library(GWASdata)
library(VariantAnnotation)
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
data(illuminaSnpADF, illuminaScanADF)
genoData <- GenotypeData(GdsGenotypeReader(gdsfile),
  scanAnnot=illuminaScanADF, snpAnnot=illuminaSnpADF)
vcf <- genoDataAsVCF(genoData, id.col="rsID")
vcf

vcffile <- tempfile()
vcfWrite(genoData, vcffile, id.col="rsID", info.cols="IntensityOnly")
vcf <- readVcf(vcffile, "hg18")
vcf
vcfCheck(genoData, vcffile, id.col="rsID")
close(genoData)
unlink(vcffile)
}
\keyword{manip}
