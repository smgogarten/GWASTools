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
\section{Changes in version 1.1.4}{
  \itemize{
    \item Added ncdfSetMissingGenotypes function.
    \item plinkCheck now writes a log file with all mismatches found.
    \item duplicateDiscordance excludes Y chrom SNPs for females.
    \item duplicateDiscordance has an option to consider only pairs
    involving the minor allele.
  }
}
\section{Changes in version 1.1.5}{
  \itemize{
    \item duplicateDiscordance has an option to compute correlation by SNP.
    \item Added scan.exclude argument to plinkCheck.
  }
}
\section{Changes in version 1.1.6}{
  \itemize{
    \item duplicateDiscordance and duplicateDiscordanceAcrossDatasets
    use only one pair of scans per subject by default.
    \item duplicateDiscordanceProbability sets small negative values to 0.
  }
}
\section{Changes in version 1.1.7}{
  \itemize{
    \item chromIntensityPlot and pseudoautoIntensityPlot have cex=0.5 by default.
    \item chromIntensityPlot colors now match anomStatsPlot colors.
    \item plinkCheck has options to skip checking parents and sex.
    \item plinkCheck sorts alleles by character to avoid phase mismatches.
    \item plinkWrite and plinkCheck print progress messages if verbose=TRUE.
  }
}
\section{Changes in version 1.1.8}{
  \itemize{
    \item Updated positions of pseudoautosomal regions.
    \item Added plinkToNcdf to convert PLINK files to NetCDF for use in GWASTools.
    \item anomSegStats checks for SNPs in centromere gaps.
  }
}
