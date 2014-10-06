\name{convertVcfGds}
\alias{convertVcfGds}
\title{Conversion from VCF to GDS}
\description{
	Extract SNP data from a VCF file
}
\usage{
convertVcfGds(vcf.filename, gds.filename, nblock=1024, compress="ZIP.max",
	      verbose=TRUE)
}
\arguments{
	\item{vcf.filename}{the file name of VCF format}
	\item{gds.filename}{the output gds file}
	\item{nblock}{the buffer lines}
	\item{compress}{the compression format for the GDS file, one of "",
	  "ZIP", "ZIP.fast", "ZIP.default", or "ZIP.max"}
	\item{verbose}{whether to show progress information}
}
\details{
	\code{convertVcfGds} extracts bi-allelic SNP genotypes from a VCF file
	and stores them in a GDS file.  All VCF rows which do not
	contain polymorphic, bi-allelic SNPs are ignored.  Unique integer IDs are generated
	for all samples and SNPs.  Sample name, SNP ID, reference and
	alternate alleles, chromosome, and position are stored in the
	GDS file as well.
	
	GDS -- Genomic Data Structures, the extended file name used for storing genetic data,
		and the file format used in the \pkg{\link{gdsfmt}} package.

	VCF -- The Variant Call Format (VCF), which is a generic format for storing DNA
		polymorphism data such as SNPs, insertions, deletions and structural variants,
		together with rich annotations.

}
\references{
	The variant call format and VCFtools.
	Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE,
	Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group.
	Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.

	\url{http://corearray.sourceforge.net/}
}
\author{Xiuwen Zheng}
\seealso{
  \code{\link{GdsGenotypeReader}}
}

\examples{
# The VCF file
vcf.file <- system.file("extdata", "sequence.vcf", package="SNPRelate")
readLines(vcf.file)

gds.file <- tempfile() 
convertVcfGds(vcf.file, gds.file)

# open GDS file
(gds <- GdsGenotypeReader(gds.file))

getScanID(gds)
getSnpID(gds)
getChromosome(gds)
getPosition(gds)
getVariable(gds, "sample.name")
getVariable(gds, "snp.rs.id")
getVariable(gds, "snp.allele")
getGenotype(gds)

# close the genotype file
close(gds)
unlink(gds.file)
}

\keyword{manip}
