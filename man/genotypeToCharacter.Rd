\name{genotypeToCharacter}
\alias{genotypeToCharacter}
\title{Convert number of A alleles to character genotypes}

\description{Converts a vector or matrix of genotypes encoded as number
  of A alleles to character strings of the form "A/B".} 

\usage{
genotypeToCharacter(geno, alleleA=NULL, alleleB=NULL, sort=TRUE)
}

\arguments{
  \item{geno}{Vector or matrix of genotype values, encoded as number of
    A alleles.  If a matrix, dimensions should be (snp, sample).}
  \item{alleleA}{Character vector with allele A.}
  \item{alleleB}{Character vector with allele B.}
  \item{sort}{Logical for whether to sort alleles lexographically ("G/T"
  instead of "T/G").}
}

\details{
  If \code{geno} is a vector, \code{alleleA} and \code{alleleB} should
  have the same length as \code{geno} or length 1 (in the latter case the values are recycled).

  If \code{geno} is a matrix, length of
  \code{alleleA} and \code{alleleB} should be equal to the number of
  rows of \code{geno}.

  If \code{alelleA} or \code{alleleB} is \code{NULL}, returned genotypes
  will have values "A/A", "A/B", or "B/B".
}

\value{
  Character vector or matrix of the same dimensions as \code{geno}.
}

\author{Stephanie Gogarten}

\seealso{\code{\link{GenotypeData}}}

\examples{
geno <- matrix(c(0,1,2,0,1,2,1,NA), nrow=4)
alleleA <- c("A","T","C","T")
alleleB <- c("C","G","G","A")
genotypeToCharacter(geno, alleleA, alleleB)
}

\keyword{manip}
