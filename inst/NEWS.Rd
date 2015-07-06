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
\section{Changes in version 1.9.4}{
  \itemize{
    \item Bug fix in getChromosome method for SnpAnnotationDataFrame
    (proper behavior of unnamed "index" argument).
  }
}
\section{Changes in version 1.9.5}{
  \itemize{
    \item Remove defunct functions.
    \item Improve efficiency of gwasExactHW, mendelErr,
    assocTestRegression (reduce number of calls to rbind).
  }
}
\section{Changes in version 1.9.6}{
  \itemize{
    \item Revised anomFilterBAF - fewer centromere spanning anomalies
    that aren't real, corrects some merging issues (previously it would
    merge sections that really were different split widths).
    Users should be aware that this will increase running time.
  }
}
\section{Changes in version 1.9.7}{
  \itemize{
    \item Added gdsSubset and gdsSubsetCheck functions to make a subset GDS file that includes only specified SNPs and scans from an existing GDS file
    \item Updated gdsImputedDosage and gdsCheckImputedDosage to account for IMPUTE2 gprobs files that have missing values (specified by three equal probability strings)
    \item Updated gdsCheckImputedDosage to produce optional logfile reporting any missing genotypes
  }
}
\section{Changes in version 1.9.8}{
  \itemize{
    \item updates/bug fixes to gdsSubset/gdsSubsetCheck - different missing value attributes may be set if sub.storage type is different.
  }
}
\section{Changes in version 1.9.9}{
  \itemize{
    \item added block size support for GDS files stored with scan,snp dimensions
    \item gdsSubset and gdsSubsetCheck now operate on the fastest dimension of the GDS file
  }
}
\section{Changes in version 1.11.1}{
  \itemize{
    \item Added ylim argument to qqPlot.
  }
}
\section{Changes in version 1.11.3}{
  \itemize{
    \item updated gdsSubset for new gdsfmt read.gdsn syntax (also changed in release version)
  }
}
\section{Changes in version 1.11.4}{
  \itemize{
    \item added optional thinThreshold argument to manhattanPlot and qqPlot functions
  }
}
\section{Changes in version 1.11.5}{
  \itemize{
    \item added pointsPerBin argument to manhattanPlot
  }
}
\section{Changes in version 1.11.6}{
  \itemize{
    \item bug fix in qqPlot, manhattanPlot when requesting thinning when bins only have 1 point
  }
}
\section{Changes in version 1.11.7}{
  \itemize{
    \item added vcfWrite function to write VCF file from GenotypeData object
  }
}
\section{Changes in version 1.11.8}{
  \itemize{
    \item added ref.allele option to vcfWrite to select either A or B as the reference allele for each SNP
  }
}
\section{Changes in version 1.11.9}{
  \itemize{
    \item changed colors of BAF plots so points can be more easily distinguished
  }
}
\section{Changes in version 1.11.10}{
  \itemize{
    \item bug fix in gwasExactHW when a block contains only 1 genotype
  }
}
\section{Changes in version 1.11.11}{
  \itemize{
    \item added vcfCheck function to compare VCF file to GenotypeData object
  }
}
\section{Changes in version 1.11.12}{
  \itemize{
    \item bug fix in assocTestRegression when a block contains only 1 SNP
  }
}
\section{Changes in version 1.11.13}{
  \itemize{
    \item in case of missing allele code, return character genotype as NA
  }
}
\section{Changes in version 1.11.14}{
  \itemize{
    \item added pasteSorted function
  }
}
\section{Changes in version 1.11.15}{
  \itemize{
    \item allowed truncate to be a numeric value or TRUE in qqPlot
  }
}
\section{Changes in version 1.11.16}{
  \itemize{
    \item expanded matching options in duplicateDiscordanceAcrossDatasets
  }
}
\section{Changes in version 1.11.17}{
  \itemize{
    \item removed duplicated .probToDosage function from ncdfImputedDosage.R source file
  }
}
\section{Changes in version 1.11.18}{
  \itemize{
    \item patch from Karl Forner to allow use of open gds objects in
    constructors for GdsReader and GdsGenotypeReader
  }
}
\section{Changes in version 1.11.19}{
  \itemize{
    \item convertNcdfGds will not write entire snp and sample
    annotations to file
    \item createDataFile replaces ncdfCreate and ncdfAddData
  }
}
\section{Changes in version 1.11.20}{
  \itemize{
    \item Fixed bug in assocTestCPH when there is no Y chromosome in the data
  }
}
\section{Changes in version 1.11.21}{
  \itemize{
    \item Fixed bug in vcfWrite to output missing data code for ID
    column
    \item Data cleaning vignette uses createDataFile instead of
    ncdfCreate and ncdfAddData
    \item Data cleaning vignette uses snpgdsOpen and snpgdsClose
  }
}
\section{Changes in version 1.11.22}{
  \itemize{
    \item convertNcdfGds and convertGdsNcdf will convert files with any
    variable names (not just genotype)
  }
}
\section{Changes in version 1.11.31}{
  \itemize{
    \item setMissingGenotypes replaces ncdfSetMissingGenotypes and gdsSetMissingGenotypes
  }
}
\section{Changes in version 1.11.32}{
  \itemize{
    \item imputedDosageFile replaces ncdfImputedDosage and gdsImputedDosage
  }
}
\section{Changes in version 1.11.33}{
  \itemize{
    \item Allow getting variables from sub-nodes in a GDS file (e.g., getVariable(GdsReader, "snp.annot/qual")).
    \item Add getNodeDescription method to GdsReader.
    \item Added examples of converting from PLINK and VCF in Formats vignette.
  }
}
\section{Changes in version 1.13.3}{
  \itemize{
    \item Bug fix in imputedDosageFile - for IMPUTE2, include columns
  from .samples file in output scan annotation.
  }
}
\section{Changes in version 1.13.6}{
  \itemize{
    \item Add getGenotypeSelection to select non-continguous SNPs and scans from GDS files.
  }
}
\section{Changes in version 1.13.7}{
  \itemize{
    \item Add dupDosageCorAcrossDatasets.
    \item Add getGenotypeSelection method to MatrixGenotypeReader.
  }
}
\section{Changes in version 1.13.8}{
  \itemize{
    \item Add drop=TRUE argument to getVariable, etc.
  }
}
\section{Changes in version 1.13.9}{
  \itemize{
    \item Add use.names argument to getGenotype and getGenotype selection
    \item Add order=c("file", "selection") argument to getGenotypeSelection
    \item duplicateDiscordanceAcrossDatasets and dupDosageCorAcrossDatasets will not match on unmapped SNPs
  }
}
\section{Changes in version 1.13.21}{
  \itemize{
    \item assocRegression replaces assocTestRegression. Only one model is allowed per function call.
    \item assocCoxPH replaces assocTestCPH. Output format is now similar to assocRegression.
    \item exactHWE replaces gwasExactHWE.
    \item assocRegression, assocCoxPH, and exactHWE include the option to select blocks of SNPs by index for easier parallelization.
    \item scan.chromosome.filter is no longer an option; use setMissingGenotypes to filter data prior to running other functions.
  }
}
\section{Changes in version 1.13.22}{
  \itemize{
    \item Removed "outfile" arguments from batchChisqTest, batchFisherTest, and mendelErr. Saving output to a file should happen outside the function calls.
    \item batchChisqTest and batchFisherTest have snp.include arguments to run on individual SNPs. Using batchFisherTest with this argument is recommended to replace the deprecated assocTestFisherExact.
  } 
}
\section{Changes in version 1.13.23}{
  \itemize{
    \item convertGdsNcdf works for transposed (sample x snp) genotype files.
  } 
}
\section{Changes in version 1.13.24}{
  \itemize{
    \item Added data tables for genome build 38: centromeres.hg38.RData, pseudoautosomal.hg38.RData, HLA.hg38.RData, pcaSnpFilters.hg38.RData.
  } 
}
\section{Changes in version 1.13.26}{
  \itemize{
    \item Added getScanAnnotation, getSnpAnnotation accessors for GenotypeData objects.
  } 
}
\section{Changes in version 1.13.27}{
  \itemize{
    \item Added "snpID" and "scanID" arguments to getGenotypeSelection.
  } 
}
\section{Changes in version 1.15.1}{
  \itemize{
    \item Added "ci" argument to qqPlot function.
  } 
}
\section{Changes in version 1.15.2}{
  \itemize{
    \item Added number of cases and controls to assocRegression output.
  } 
}
\section{Changes in version 1.15.3}{
  \itemize{
    \item Added beta and standard error for GxE term to assocRegression output.
  } 
}
\section{Changes in version 1.15.5}{
  \itemize{
    \item Added option to read genotypes coded with nucleotides in createDataFile.
  } 
}
\section{Changes in version 1.15.7}{
  \itemize{
    \item Added option to reorder samples in vcfWrite.
  }   
}
\section{Changes in version 1.15.8}{
  \itemize{
    \item alleleFrequency includes scans with missing sex.
  }   
}
\section{Changes in version 1.15.10}{
  \itemize{
    \item createDataFile converts non-finite values to NA.
  }   
}
\section{Changes in version 1.15.11}{
  \itemize{
    \item added col argument to manhattanPlot
  }   
}
\section{Changes in version 1.15.12}{
  \itemize{
    \item allow character scanID in createDataFile
  }   
}
