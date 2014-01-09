
# SnpAnnotationDataFrame
# Snp Annotation stored as an AnnotatedDataFrame
setClass("SnpAnnotationDataFrame",
         contains = "AnnotatedDataFrame",
         representation(idCol = "character",
                        chromosomeCol = "character",
                        positionCol = "character",
                        alleleACol = "character",
                        alleleBCol = "character",
                        autosomeCode = "integer",
                        XchromCode = "integer",
                        YchromCode = "integer",
                        XYchromCode = "integer",
                        MchromCode = "integer"),
         prototype(idCol = "snpID",
                   chromosomeCol = "chromosome",
                   positionCol = "position",
                   alleleACol = "alleleA",
                   alleleBCol = "alleleB",
                   autosomeCode = 1:22,
                   XchromCode = 23L,
                   YchromCode = 25L,
                   XYchromCode = 24L,
                   MchromCode = 26L))

# SnpAnnotationSQLite
# Snp Annotation stored in a SQLite database
setClass("SnpAnnotationSQLite",
         representation(dbpath = "character",
                        dbcon = "SQLiteConnection",
                        annotationTable = "character",
                        metadataTable = "character",
                        idCol = "character",
                        chromosomeCol = "character",
                        positionCol = "character",
                        alleleACol = "character",
                        alleleBCol = "character",
                        autosomeCode = "integer",
                        XchromCode = "integer",
                        YchromCode = "integer",
                        XYchromCode = "integer",
                        MchromCode = "integer"),
         prototype(idCol = "snpID",
                   chromosomeCol = "chromosome",
                   positionCol = "position",
                   alleleACol = "alleleA",
                   alleleBCol = "alleleB",
                   autosomeCode = 1:22,
                   XchromCode = 23L,
                   YchromCode = 25L,
                   XYchromCode = 24L,
                   MchromCode = 26L))

# SnpAnnotationReader
# Generic reader class
# Add more Reader classes as necessary
# NULL is included so GenotypeData, etc. are not required to have annotation
setClassUnion("SnpAnnotationReader", c("SnpAnnotationDataFrame",
                                       "SnpAnnotationSQLite", "NULL"))


# ScanAnnotationDataFrame
# Scan Annotation stored as an AnnotatedDataFrame
setClass("ScanAnnotationDataFrame",
         contains = "AnnotatedDataFrame",
         representation(idCol = "character",
                        sexCol = "character"),
         prototype(idCol = "scanID",
                   sexCol = "sex"))

# ScanAnnotationSQLite
# Scan Annotation stored in a SQLite database
setClass("ScanAnnotationSQLite",
         representation(dbpath = "character",
                        dbcon = "SQLiteConnection",
                        annotationTable = "character",
                        metadataTable = "character",
                        idCol = "character",
                        sexCol = "character"),
         prototype(idCol = "scanID",
                   sexCol = "sex"))

# ScanAnnotationReader
# Generic reader class
# Add more Reader classes as necessary
# NULL is included so GenotypeData, etc. are not required to have annotation
setClassUnion("ScanAnnotationReader", c("ScanAnnotationDataFrame",
                                        "ScanAnnotationSQLite", "NULL"))


# NcdfReader
# Wrapper for the ncdf library - holds ncdf file handler
# ncdf is S3 class, so create S4 virtualization
setOldClass("ncdf")
setClass("NcdfReader",
         representation(filename = "character",
                        handler = "ncdf"))

# NcdfGenotypeReader
# Reads genotype data stored in netCDF format
setClass("NcdfGenotypeReader",
         contains = "NcdfReader",
         representation(snpDim = "character",
                        scanDim = "character",
                        snpIDvar = "character",
                        chromosomeVar = "character",
                        positionVar = "character",
                        scanIDvar = "character",
                        genotypeVar = "character",
                        autosomeCode = "integer",
                        XchromCode = "integer",
                        YchromCode = "integer",
                        XYchromCode = "integer",
                        MchromCode = "integer"),
         prototype(snpDim = "snp",
                   scanDim = "sample",
                   snpIDvar = "snp",
                   chromosomeVar = "chromosome",
                   positionVar = "position",
                   scanIDvar = "sampleID",
                   genotypeVar = "genotype",
                   autosomeCode = 1:22,
                   XchromCode = 23L,
                   YchromCode = 25L,
                   XYchromCode = 24L,
                   MchromCode = 26L))

# MatrixGenotypeReader
# Stores genotype data in a matrix
setClass("MatrixGenotypeReader",
         representation(snpID = "integer",
                        chromosome = "integer",
                        position = "integer",
                        scanID = "integer",
                        genotype = "matrix",
                        autosomeCode = "integer",
                        XchromCode = "integer",
                        YchromCode = "integer",
                        XYchromCode = "integer",
                        MchromCode = "integer"),
         prototype(autosomeCode = 1:22,
                   XchromCode = 23L,
                   YchromCode = 25L,
                   XYchromCode = 24L,
                   MchromCode = 26L))

# GdsReader
# holds GDS file handler (gdsfmt library)
setOldClass("gds.class")
setClass("GdsReader",
         representation(filename = "character",
                        handler = "gds.class"))

# GdsGenotypeReader
# Reads genotype data stored in GDS format
setClass("GdsGenotypeReader",
         contains = "GdsReader",
         representation(snpIDvar = "character",
                        chromosomeVar = "character",
                        positionVar = "character",
                        scanIDvar = "character",
                        genotypeVar = "character",
                        alleleVar = "character",
                        autosomeCode = "integer",
                        XchromCode = "integer",
                        YchromCode = "integer",
                        XYchromCode = "integer",
                        MchromCode = "integer",
                        genotypeDim = "character",
                        missingValue = "integer"),
         prototype(snpIDvar = "snp.id",
                   chromosomeVar = "snp.chromosome",
                   positionVar = "snp.position",
                   scanIDvar = "sample.id",
                   genotypeVar = "genotype",
                   alleleVar = "snp.allele",
                   autosomeCode = 1:22,
                   XchromCode = 23L,
                   YchromCode = 25L,
                   XYchromCode = 24L,
                   MchromCode = 26L,
                   genotypeDim = "",
                   missingValue = 3))

# GenotypeReader
# Generic reader class
# Add more Reader classes as necessary
setClassUnion("GenotypeReader", c("NcdfGenotypeReader", "MatrixGenotypeReader",
                                  "GdsGenotypeReader"))

# GenotypeData
setClass("GenotypeData",
         representation(data = "GenotypeReader",
                        snpAnnot = "SnpAnnotationReader",
                        scanAnnot = "ScanAnnotationReader"),
         prototype(snpAnnot = NULL,
                   scanAnnot = NULL))

# NcdfIntensityReader
# Reads intensity data stored in netCDF format
setClass("NcdfIntensityReader",
         contains = "NcdfReader",
         representation(snpDim = "character",
                        scanDim = "character",
                        snpIDvar = "character",
                        chromosomeVar = "character",
                        positionVar = "character",
                        scanIDvar = "character",
                        qualityVar = "character",
                        xVar = "character",
                        yVar = "character",
                        bafVar = "character",
                        lrrVar = "character",
                        autosomeCode = "integer",
                        XchromCode = "integer",
                        YchromCode = "integer",
                        XYchromCode = "integer",
                        MchromCode = "integer"),
         prototype(snpDim = "snp",
                   scanDim = "sample",
                   snpIDvar = "snp",
                   chromosomeVar = "chromosome",
                   positionVar = "position",
                   scanIDvar = "sampleID",
                   qualityVar = "quality",
                   xVar = "X",
                   yVar = "Y",
                   bafVar = "BAlleleFreq",
                   lrrVar = "LogRRatio",
                   autosomeCode = 1:22,
                   XchromCode = 23L,
                   YchromCode = 25L,
                   XYchromCode = 24L,
                   MchromCode = 26L))

# IntensityReader
# Generic reader class
# Add more Reader classes as necessary
setClassUnion("IntensityReader", c("NcdfIntensityReader"))

# IntensityData
setClass("IntensityData",
         representation(data = "IntensityReader",
                        snpAnnot = "SnpAnnotationReader",
                        scanAnnot = "ScanAnnotationReader"),
         prototype(snpAnnot = NULL,
                   scanAnnot = NULL))
