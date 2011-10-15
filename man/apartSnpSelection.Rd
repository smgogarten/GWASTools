\name{apartSnpSelection}
\alias{apartSnpSelection}
\title{Random selection of SNPs
}
\description{Randomly selects SNPs for which each pair is at least as far apart as the specified basepair distance.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usage{
apartSnpSelection(chromosome, position, min.dist = 1e+05, 
                  init.sel = NULL,  max.n.chromosomes = -1, 
                  verbose = TRUE)
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{chromosome}{An integer vector containing the chromosome for each
  SNP.  Valid values are 1-26, any other value will be interpreted as
  missing and not selected.
}
  \item{position}{A numeric vector of the positions (in basepairs) of the SNPs.
}
  \item{min.dist}{A numeric value to specify minimum distance required (in basepairs).
}
  \item{init.sel}{A logical vector indicating the initial SNPs to be included.
}
  \item{max.n.chromosomes}{A numeric value specifying the maximum number of SNPs to return per chromosome, "-1" means no number limit.
}
  \item{verbose}{A logical value specifying whether to show progress information while running.
}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\details{
  \code{apartSnpSelection} selects SNPs randomly with the condition that
  they are at least as far apart as \code{min.dist} in basepairs.  The
  starting set of SNPs can be specified with \code{init.sel}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\value{
A logical vector indicating which SNPs were selected.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\author{Xiuwen Zheng
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
library(GWASdata)
data(affy_snp_annot)
pool <- affy_snp_annot$chromosome < 23
rsnp <- apartSnpSelection(affy_snp_annot$chromosome, affy_snp_annot$position,
                          min.dist=15000, init.sel=pool)
}
\keyword{manip}

