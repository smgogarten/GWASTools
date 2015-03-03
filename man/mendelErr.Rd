\name{mendelErr}
\alias{mendelErr}
\title{Mendelian Error Checking}
\description{Mendelian and mtDNA inheritance tests.}

\usage{
mendelErr(genoData, mendel.list, snp.exclude=NULL,
          error.by.snp=TRUE, error.by.snp.trio=FALSE,
          verbose=TRUE)
}
\arguments{
	\item{genoData}{\code{\link{GenotypeData}} object, must have
	  scan variable "sex"}
	\item{mendel.list}{A \code{\link{mendelList}} object, to specify trios.}
	\item{snp.exclude}{An integer vector with snpIDs of SNPs to
	  exclude.  If \code{NULL} (default), all SNPs are used.}
	\item{error.by.snp}{Whether or not to output Mendelian errors per SNP. This will only return the total number of trios checked
and the total number of errors for each SNP. The default value is \code{TRUE}.}
	\item{error.by.snp.trio}{Whether or not to output Mendelian errors per SNP for each trio. This will return the total number of trios checked
and the total number of errors for each SNP as well as indicators of which SNPs have an error for each trio.
The default value is \code{FALSE}.
NOTE:  \code{error.by.snp} must be set to \code{TRUE} as well in order to use this option.
NOTE:  Using this option causes the output to be very large that may be slow to load into R.
}
	\item{verbose}{If \code{TRUE} (default), will print status updates while the function runs.}
}

\details{
\code{genoData} must contain the scan annotation variable "sex".
Chromosome index: 1..22 autosomes, 23 X, 24 XY, 25 Y, 26 mtDNA, 27
missing.
}

\value{
\code{mendelErr} returns an object of class "mendelClass".
The object contains two data
frames: "trios" and "all.trios", and a list: "snp" (if \code{error.by.snp} is specified to
be \code{TRUE}). If there are no duplicate samples in the dataset, "trios" will be the same as
"all.trios". Otherwise, "all.trios" contains the results of all combinations of
duplicate samples, and "trios" only stores the average values of unique trios. i.e: "trios" averages duplicate
samples for each unique subject trio.
"trios" and "all.trios" contain the following components:

	\item{fam.id}{Specifying the family ID from the mendel.list object used as input.}
	\item{child.id}{Specifying the offspring ID from the mendel.list object used as input.}
	\item{child.scanID}{Specifying the offspring scanID from the mendel.list object used as input. (only in "all.trios")}
	\item{father.scanID}{Specifying the father scanID from the mendel.list object used as input. (only in "all.trios")}
	\item{mother.scanID}{Specifying the mother scanID from the mendel.list object used as input. (only in "all.trios")}
	\item{Men.err.cnt}{The number of SNPs with Mendelian errors in this trio.}
	\item{Men.cnt}{The total number of SNPs checked for Mendelian errors in this trio.
		It excludes those cases where the SNP is missing in the offspring and
		those cases where it is missing in both parents.
		Hence, Mendelian error rate = \code{Men.err.cnt} / \code{Men.cnt}.}
	\item{mtDNA.err}{The number of SNPs with mtDNA inheritance errors in this trio.}
	\item{mtDNA.cnt}{The total number of SNPs checked for mtDNA inheritance errors in this trio.
		It excludes those cases where the SNP is missing in the offspring
		and in the mother. Hence, mtDNA error rate = \code{mtDNA.err} / \code{mtDNA.cnt} .}
	\item{chr1, ..., chr25}{The number of Mendelian errors in each chromosome for this trio.}
	
"snp" is a list that contains the following components:	

	\item{check.cnt}{A vector of integers, indicating the number of trios valid for checking on each SNP.}
	\item{error.cnt}{A vector of integers, indicating the number of trios with errors on each SNP.}
	\item{familyid.childid}{A vector of indicators (0/1) for whether or not any of the duplicate trios for the unique trio,
	"familyid.childid", have a Mendelian error on each SNP. (Only if \code{error.by.snp.trio} is specified to be \code{TRUE}).}
}

\author{Xiuwen Zheng, Matthew P. Conomos}

\seealso{\code{\link{mendelList}}}

\examples{
library(GWASdata)
data(illuminaScanADF)
scanAnnot <- illuminaScanADF

# generate trio list
men.list <- mendelList(scanAnnot$family, scanAnnot$subjectID,
  scanAnnot$father, scanAnnot$mother, scanAnnot$sex,
  scanAnnot$scanID)

# create genoData object
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(gdsfile)
genoData <- GenotypeData(gds, scanAnnot=scanAnnot)

# Run!
R <- mendelErr(genoData, men.list, error.by.snp.trio = TRUE)

names(R)
# [1] "trios"     "all.trios" "snp"

names(R$trios)
#  [1] "fam.id"      "child.id"    "Men.err.cnt" "Men.cnt"     "mtDNA.err"  
#  [6] "mtDNA.cnt"   "chr1"        "chr2"        "chr3"        "chr4"       
# [11] "chr5"        "chr6"        "chr7"        "chr8"        "chr9"       
# [16] "chr10"       "chr11"       "chr12"       "chr13"       "chr14"      
# [21] "chr15"       "chr16"       "chr17"       "chr18"       "chr19"      
# [26] "chr20"       "chr21"       "chr22"       "chr23"       "chr24"      
# [31] "chr25"

# Mendelian error rate = Men.err.cnt / Men.cnt
data.frame(fam.id = R$trios$fam.id, child.id = R$trios$child.id,
           Mendel.err.rate = R$trios$Men.err.cnt / R$trios$Men.cnt)

names(R$snp)
summary(R$snp$check.cnt)

# summary Mendelian error for first family
summary(R$snp[[1]])

close(genoData)
}

\keyword{Mendelian}

