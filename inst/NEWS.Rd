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
  }
}
\section{Changes in version 1.1.9}{
  \itemize{
    \item anomSegStats checks for SNPs in centromere gaps.
    \item anomStatsPlot has option to plot LRR/BAF individually (for
    greater flexibility in layout).
    \item Updates to arguments for plot titles in chromIntensityPlot,
    anomStatsPlot, and pseudoautoIntensityPlot for consistency.
    \item plinkCheck has map.alt argument to override default
    GenotypeData -> PLINK annotation conversion.
  }
}
\section{Changes in version 1.3.1}{
  \itemize{
    \item close method of NcdfReader returns invisibly.
  }
}
\section{Changes in version 1.3.2}{
  \itemize{
    \item Use lazy loading of data.
    \item manhattanPlot and snpCorrelationPlot accept character vectors
    of chromosome; chrom.labels argument no longer used.
  }
}
\section{Changes in version 1.3.3}{
  \itemize{
    \item Updated vignette.
  }
}
\section{Changes in version 1.3.4}{
  \itemize{
    \item Added chromosome ideograms to chromIntensityPlot and
    anomStatsPlot.  anomStatsPlot has an option to put multiple
    anomalies on the same plot.    
  }
}
\section{Changes in version 1.3.5}{
  \itemize{
    \item Added function asSnpMatrix to convert a GenotypeData object
    to a SnpMatrix object for use with snpStats.
  }
}
\section{Changes in version 1.3.6}{
  \itemize{
    \item Removed automatic subtitle from qqPlot.
    \item Allow selection of theoretical boundaries to draw in ibdPlot.
  }
}
\section{Changes in version 1.3.7}{
  \itemize{
    \item In batchChisqTest, Yates correction cannot be bigger than the
    terms it corrects.  Changed to match bug fix to chisq.test in R 2.15.1.
  }
}
\section{Changes in version 1.3.8}{
  \itemize{
    \item Pedigree functions accept non-numeric identifiers and provide
    additional output.
  }
}
\section{Changes in version 1.3.9}{
  \itemize{
    \item Fixed bug in assocTestCPH that could lead to false positives
    if additive model failed but GxE model did not.
    \item Allow multiple variables for stratified analysis in assocTestCPH.
  }
}
