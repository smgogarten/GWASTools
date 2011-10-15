\name{relationsMeanVar}
\alias{relationsMeanVar}
\docType{data}
\title{Mean and Variance information for full-sibs, half-sibs, first-cousins}
\description{
  Computes theoretical mean and covariance matrix for k0 vs. k1 ibd coefficients
  for full-sib relationship along with inverse and eigenvalues/vectors of
  the covariance matrix. 

Computes theoretical means and variances for half-sib 
  relationship and for first-cousin relationship. 
}
\usage{relationsMeanVar}
\format{
  A list with the following entries:
 \describe{
     \item{\code{FullSibs}}{list with following entries:
       \itemize{
          \item \code{mean}:  mean of (k0,k1) for full-sibs
          \item \code{cov}:  covariance matrix for full-sibs
          \item \code{invCov}:  inverse of the covariance matrix
         \item \code{eigvals}:  eigenvalues of the inverse covariance matrix
         \item \code{eigvectors}:  eigenvectors of the inverse covariance matrix
       }
     }
     \item{\code{HalfSibs}}{list with following entries:
         \itemize{
            \item \code{mean}:  mean of (k0,k1) for half-sibs
            \item \code{var}:  variance for half-sibs
         }
      }
     \item{\code{FirstCousins}}{list with following entries:
         \itemize{
            \item \code{mean}:  mean of (k0,k1) for first-cousins
            \item \code{var}:  variance for first-cousin
         }
      } 
  }    
}
\source{computed by Cecelia Laurie using the referenced papers}
\references{
  Hill, W.G. and B.S. Weir (2011) Variation in actual relationship as a consequence 
 of Mendelian sampling and linkage,\emph{Genet. Res., Camb.}, \bold{93}, 47--64.

Kong, X., \emph{et al} (2004) A combined physical-linkage map of the human genome,
\emph{American Journal of Human Genetics}, \bold{75}, 1143--1148.

}
\examples{
data(relationsMeanVar)
FS<-relationsMeanVar$FullSibs
FScov<-FS$cov   #gives covariance matrix for full-sibs
HS<-relationsMeanVar$HalfSibs
HSvar<-HS$var  #gives variance for half-sibs
}
\keyword{datasets}
