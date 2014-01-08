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
\section{Changes in version 1.3.10}{
  \itemize{
    \item Fixed bug in chromIntensityPlot so ideogram scales correctly
    if SNPs are excluded.
  }
}
\section{Changes in version 1.3.11}{
  \itemize{
    \item Fixed bug in duplicateDiscordance when Y chromosome is not
    included.
  }
}
\section{Changes in version 1.3.12}{
  \itemize{
    \item Added GdsReader and GdsGenotypeReader classes with dependency
    on gdsfmt.  GenotypeData objects can also be created with
    GdsGenotypeReader objects in the "data" slot.
  }
}
\section{Changes in version 1.3.13}{
  \itemize{
    \item Added vignette on converting VCF to NetCDF with annotation.
    \item Prevent duplicateDiscordance from checking correlation by SNP
    in cases of no variation.
  }
}
\section{Changes in version 1.3.14}{
  \itemize{
    \item assocTestRegression returns NA for snps where cases or controls are monomorphic, added assocTestFisherExact to use in that case.
    \item Added snp.exclude argument to pseudoautoIntensityPlot.
    \item Bug fix in messages reporting file read times when creating or
    checking netCDF files.
  }
}
\section{Changes in version 1.3.15}{
  \itemize{
    \item Bug fix in pedigreePairwiseRelatedness related to use of character identifiers.
  }
}
\section{Changes in version 1.3.16}{
  \itemize{
    \item Added convertVcfGds to extract bi-allelic SNPs from a VCF
    file.
    \item Added ncdfImputedDosage to convert output from common
    imputation programs to NetCDF.  assocTestRegression has an
    additional argument dosage=TRUE to be used with these files.
    \item Added vignette describing GWASTools data structures.
  }
}
\section{Changes in version 1.5.1}{
  \itemize{
    \item Added option to duplicateDiscordanceAcrossDatasets to
    calculate minor allele discordance.
  }
}
\section{Changes in version 1.5.2}{
  \itemize{
    \item Added "getAlleleA" and "getAlleleB" methods to SnpAnnotation*
    and GenotypeData objects.
    \item Added genotypeToCharacter function to convert genotypes from
    number of A alleles to A/B format.
    \item getGenotype for GenotypeData has option char=TRUE to return
    character genotypes in A/B format.
  }
}
\section{Changes in version 1.5.3}{
  \itemize{
    \item Removed "alleleA.col" and "alleleB.col" options from plink
    functions, as "alleleA" and "alleleB" are now standard names.
    \item Added "getAlleleA" and "getAlleleB" methods to
    GdsGenotypeReader.
    \item Added "getDimension" method to NcdfReader.
  }
}
\section{Changes in version 1.5.4}{
  \itemize{
    \item Removed hard-coding of autosomes as 1:22; can now set a vector
    of integer codes corresponding to autosomes
    with "autosomeCode" argument at object creation and retrieve with
    "autosomeCode" methods.  This change makes GWASTools compatible with
    non-human organisms.
    \item Added option to duplicateDiscordanceAcrossDatasets to
    count missing data as discordance.
    \item Added option to start axes of genoClusterPlot at 0.
  }
}
\section{Changes in version 1.5.5}{
  \itemize{
    \item Deprecated pedigreeClean and pedigreeFindDuplicates.
    pedigreeCheck now encompasses all pedigree checks and should be used
    instead.
    \item Added pedigreeMaxUnrelated to find the maximum set of
    unrelated members of a pedigree.
    \item Added additional output column "MAF" to matrix returned by
    alleleFrequency.
  }
}
\section{Changes in version 1.5.6}{
  \itemize{
    \item Added function minorAlleleSensitivitySpecificity.
  }
}
\section{Changes in version 1.5.7}{
  \itemize{
    \item Fixed bug in genotypeToCharacter that resulted in calls to
    getGenotype(char=TRUE) for a single SNP to return NA.
    \item Renamed minorAlleleSensitivitySpecificity to
    minorAlleleDetectionAccuracy and added additional output.
  }
}
\section{Changes in version 1.5.8}{
  \itemize{
    \item Adding missing value support to GdsReader.
    \item Fixed bug in getAttribute method for GdsReader.
    \item Updated GdsReader for compatibility with gdsfmt 0.9.11 (no
    longer compatible with older versions).
  }
}
\section{Changes in version 1.5.9}{
  \itemize{
    \item assocTestRegression computes allele counts separately for each
    model.
    \item convertNcdfGds uses information from a SnpAnnotationDataFrame to store allele
    and chromosome codes in the GDS file.
  }
}
\section{Changes in version 1.7.3}{
  \itemize{
    \item Changed labeling of IBD plots from "HS" to "Deg2" and "FC" to "Deg3."
    \item Bug fix in pedigreePairwiseRelatedness - no more warning about
    multiple values passed to stringsAsFactor.
    \item pedigreeClean and pedigreeFindDuplicates are now defunct.  Use
    pedigreeCheck instead.
  }
}
\section{Changes in version 1.7.4}{
  \itemize{
    \item More informative error messages in anomDetectBAF and
    anomDetectLOH.
  }
}
\section{Changes in version 1.7.5}{
  \itemize{
    \item gwasExactHW will run on all chromosomes except (Y,M), rather
    than (autosome,X,XY) only.
  }
}
\section{Changes in version 1.7.6}{
  \itemize{
    \item Updated vignette for SNPRelate 0.9.16.
  }
}
\section{Changes in version 1.7.7}{
  \itemize{
    \item Changed handling of GxE interaction variables in assocTestRegression.
  }
}
\section{Changes in version 1.7.8}{
  \itemize{
    \item Added gdsSetMissingGenotypes, updated argument names in ncdfSetMissingGenotypes.
    \item Changed colorscheme in manhattanPlot.R.
    \item Bug fix in ibdPlot - diagonal orange bars are back.
    \item Bug fix in plinkWrite for writing just one sample.
    \item Bug fix in printing pedigreeCheck error message.
  }
}
\section{Changes in version 1.9.1}{
  \itemize{
    \item Documentation now located in vignettes/ folder.
    \item Added ibdAssignRelationshipsKing.
    \item Added support for genotype GDS files with scan x snp dimensions in GdsGenotypeReader.
  }
}
\section{Changes in version 1.9.2}{
  \itemize{
    \item ScanAnnotationDataFrame and ScanAnnotationSQLite allow
    non-integer scanID.
    \item Fix getAlleleA and getAlleleB in GdsGenotypeReader to work
    with indels.
  }
}
\section{Changes in version 1.9.3}{
  \itemize{
    \item Added gdsImputedDosage function.
    \item GdsGenotypeReader can return transposed genotypes.
  }
}

