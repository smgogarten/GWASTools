\name{duplicateDiscordanceProbability}
\alias{duplicateDiscordanceProbability}
\title{Probability of duplicate discordance}
\description{
\code{duplicateDiscordanceProbability} calculates the probability of
observing discordant genotypes for duplicate samples.
}

\usage{
duplicateDiscordanceProbability(npair,
                                error.rate = c(1e-5, 1e-4, 1e-3, 1e-2),
                                max.disc = 7)
}

\arguments{
  \item{npair}{The number of pairs of duplicate samples.}
  \item{error.rate}{A numeric vector of error rates (i.e., the rate at
    which a genotype will be called incorrectly).}
  \item{max.disc}{The maximum number of discordances for which to
  compute the probability.}
}

\details{
  Since there are three possible genotypes, one call is correct and the
  other two are erroneous, so theoretically there are two error rates, a
  and b. The probability that duplicate genotyping instances of the same
  subject will give a discordant genotype is 2[(1 - a - b)(a + b) +
  ab]. When a and b are very small, this is approximately 2(a + b) or
  twice the total error rate.
  This function assumes that a == b, and the argument \code{error.rate}
  is the total error rate a + b.

  Any negative values for the probability (due to precision problems for
  very small numbers) are set to 0.
}

\value{
  This function returns a matrix of probabilities, where the column
  names are error rates and the row names are expected number of
  discordant genotypes (>0 through >\code{max.disc}).
}

\author{Cathy Laurie}

\seealso{\code{\link{duplicateDiscordance}},
  \code{\link{duplicateDiscordanceAcrossDatasets}} }
  
\examples{
disc <- duplicateDiscordanceProbability(npair=10, error.rate=c(1e-6, 1e-4))

#probability of observing >0 discordant genotypes given an error rate 1e-6
disc[1,1]

#probability of observing >1 discordant genotypes given an error rate 1e-4
disc[2,2] 
}

\keyword{distributiion}
