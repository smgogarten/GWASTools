\name{pedigreeMaxUnrelated}
\alias{pedigreeMaxUnrelated}

\title{
 Find a maximal set of unrelated individuals in a subset of a pedigree.
}
\description{
Given a full pedigree (with no duplicates and no one-person families), this function finds a maximal set of unrelated individuals in a specified subset of the pedigree.
This is done family by family. The full pedigree is checked for inconsistencies and an error message is given if inconsistencies are found (see \code{\link{pedigreeCheck}}).
Maximal sets are not unique; there is an option for the user to identify preference(s) in the choice of individuals.
}
\usage{
pedigreeMaxUnrelated(pedigree, pref = NULL)
}

\arguments{
  \item{pedigree}{
A dataframe containing the full pedigree with columns 'family', 'individ', 'mother', 'father', 'sex', and 'selset'.
The variables 'family', 'individ', 'mother', 'father' contain the identifiers for family, individual, individual's mother and individual's father.
Identifiers can be integer, numeric or character but identifiers for mother and father for founders are assumed to be 0. 
The variable 'sex' contains the individual's sex (coded as "M" or "F"). 
The varible 'selset' is coded as 1 = if individual is in the subset of interest and 0 otherwise.
The dataframe can contain an optional variable indicating preferences for choosing individuals. See the item \code{pref} below. 
}
  \item{pref}{
\code{pref} = the name of the (optional) preference column in \code{samp}. Preferences can be layered.
This variable must have integer  or numeric values greater than or equal to 1 where a lower value indicates higher preference. 
If \code{pref} is missing, the default is to prefer choosing founders. 
}
}
\details{
Commonly used for selecting a maximal unrelated set of genotyped individuals from a pedigree ('selset' = 1 if individual is genotyped and 0 otherwise).

An example of the use of a layered preference variable: if one wanted to prefer cases over controls and then prefer founders, 
the preference variable would = 1 for cases, 2 = founder, 3 = otherwise.
}
\value{
A dataframe with variables 'family' = family identifier and 'Individ' = individual identifier of individuals in the maximal unrelated set.
}

\author{
Cecelia Laurie
}
\note{
Since \code{pedigreeMaxUnrelated} does not accept one-person families included in the input pedigree,
to get a complete maximal set of unrelated individuals from a specified subset of the pedigree, the user
will need to append to the output from the function the one-person family (singleton) individuals from the specified subset.
}


\seealso{
\code{\link{pedigreeCheck}}, \code{\link{pedigreePairwiseRelatedness}}
}
\examples{
## Example set 1
family <- rep("A",8)
individ <- c("a","b","c","d","e","f","g","h")
mother <- c(0,"a","b",0,"f",0,0,"f")
father <- c(0,"d","e",0,"g",0,0,"g")
sex <- c(rep("F",3),"M","M","F","M","F")
pedigree <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)

## preference default (i.e. choose founders if possible)
pedigree$selset <- 1  # all selected
pedigreeMaxUnrelated(pedigree)   # chose the founders
#  family Individ
#1      A       a
#2      A       d
#3      A       f
#4      A       g

sel <- is.element(pedigree$individ,c("a","f","g"))
pedigree$selset[sel] <- 0  #only one founder 'd' in desired subset

#  default preference of founders
pedigreeMaxUnrelated(pedigree)
#  family Individ
#1      A       d    #founder
#2      A       e

## preference choice
pedigree$pref <- 2
sel2 <- is.element(pedigree$individ, c("c","h")) # preferred choices 
pedigree$pref[sel2] <- 1
pedigreeMaxUnrelated(pedigree,pref="pref")
#  family Individ
#1      A       h
#2      A       b

## add preference layer of secondary choice of founders
pedigree$pref <- 3
sel2 <- pedigree$mother==0 & pedigree$father==0
sel1 <- is.element(pedigree$individ, c("c","h"))
pedigree$pref[sel2] <- 2
pedigree$pref[sel1] <- 1
pedigreeMaxUnrelated(pedigree,pref="pref") 
#  family Individ
#1      A       h   #top pref
#2      A       d   #founder
#Note that the other top preference 'c' is related to everyone so not chosen

## Example Set 2
family <- c(1,1,1,1,2,2,2,2,2)
individ <- c(2,1,3,4,"A5","A6","A7","A8","A9")
mother <- c(3,3,0,0,0,0,"A5","A5",0)
father <- c(4,4,0,0,0,0,"A6","A9",0)
sex <- c("F","M","F","M","F","M","M","M","M")
pedigree <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
pedigree$selset <- 1
pedigree$selset[is.element(pedigree$individ, c("A5",4))] <- 0
pedigree$pref <- 2
pedigree$pref[is.element(pedigree$individ,c("A8","A7"))] <- 1
pedigreeMaxUnrelated(pedigree,pref="pref") 
#  family Individ
#1      1       2
#2      2      A6
#3      2      A8
# NOTE: in using the pref option there is NO preference for family 1  
# so will select one unrelated from family 1: 
# individual 2 is selected since it is first in selset to be listed in pedigree   

pedigree$pref <- 2
pedigree$pref[is.element(pedigree$individ,c("A8","A7"))] <- 1
sel <- pedigree$family==1 & pedigree$mother==0 & pedigree$father==0 #founders
pedigree$pref[sel] <- 1
pedigreeMaxUnrelated(pedigree,pref="pref") 
#  family Individ
#1      1       3
#2      2      A6
#3      2      A8
}

\keyword{manip }

