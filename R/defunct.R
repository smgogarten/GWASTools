assocTestRegression <- function(genoData,
                                outcome,
                                model.type,
                                covar.list = NULL,
                                ivar.list = NULL,
                                gene.action.list = NULL,
                                dosage = FALSE,
                                scan.chromosome.filter = NULL,
                                scan.exclude = NULL,
                                CI = 0.95,
                                robust = FALSE,
                                LRtest = TRUE,
                                chromosome.set = NULL,
                                block.set = NULL,
                                block.size = 5000,
                                verbose = TRUE,
                                outfile = NULL){

    .Defunct("assocRegression")
}


assocTestFisherExact <- function(dat, outfile = NULL){

    .Defunct("batchFisherTest")
}


assocTestCPH <- function(
	genoData,	# GenotypeData object containing sex and phenotypes
	event,	# name of variable in genoData for event to analyze
	time.to.event,		# name of variable in genoData for time to event
	covars,		# vector of covariate terms for model (can include interactions as 'a:b', main effects correspond to variable names in genoData)
	factor.covars=NULL,		# vector of names of covariates to be converted to factor
	scan.chromosome.filter = NULL,  # matrix of T/F for scan by chromosome for chromosome-specific selection of samples
        scan.exclude = NULL,
	maf.filter = FALSE,  # whether to filter results returned using maf  > 75/2n where n = number of events
	GxE = NULL,     # name of the covariate to use for E if genotype-by-environment (i.e. SNP:E) model is to be analyzed, in addition to the main effects (E can be a covariate interaction)
	strata.vars = NULL,  # name of variable to stratify on for a stratified analysis (use NULL if no stratified analysis needed)
	chromosome.set = NULL, 	# vector of chromosome numbers (corresponding to format of "chromosome" in genoData - i.e. integer codes)
	block.size = 5000,	# number of SNPs from a given chromosome to read in one block from genoData
        verbose = TRUE,
        outfile = NULL
){

    .Defunct("assocCoxPH")
}


gwasExactHW <- function (genoData,
                         scan.chromosome.filter = NULL,
                         scan.exclude = NULL,
                         geno.counts = TRUE,
                         chromosome.set = NULL,
                         block.size = 5000,                      
                         verbose = TRUE,
                         outfile = NULL) 
{
    .Defunct("exactHWE")
}
