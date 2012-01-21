\name{NEWS}
\title{NEWS for GWASTools}
\section{Changes in version 1.1.1}{
  \itemize{
    \item Added functions plinkWrite and plinkCheck for writing and
    checking PLINK ped and map files.
    \item Added pcaSnpFilters data set for identifying regions with high
    PC-SNP correlation.
  }
}
\section{Changes in version 1.1.2}{
  \itemize{
    \item Added LR tests to assocTestRegression.
    \item Bug fix in calculation of mean odds ratio in batchFisherTest.
    \item Bug fix in missingGenotypeByScanChrom for data sets with only
    one female.
  }
}
\section{Changes in version 1.1.3}{
  \itemize{
    \item batchChisqTest and batchFisherTest now return n results for n
    batches even if n=2.
    \item batchFisherTest has return.by.snp=FALSE as default.
  }
}
