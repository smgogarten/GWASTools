test_makeFilter <- function() {
    snp <- data.frame(snpID=1:3,
                      chromosome=1:3,
                      position=1:3,
                      filt1=c(TRUE, TRUE, FALSE),
                      filt2=c(TRUE, FALSE, FALSE))
    snpAnnot <- SnpAnnotationDataFrame(snp)
    varMetadata(snpAnnot)[c("filt1", "filt2"), "labelDescription"] <-
        c("filter 1", "filter 2")
    mgr <- MatrixGenotypeReader(matrix(0,nrow=3,ncol=2),
                                snpID=snp$snpID,
                                chromosome=snp$chromosome,
                                position=snp$position,
                                scanID=1:2)
    genoData <- GenotypeData(mgr, snpAnnot=snpAnnot)

    f <- GWASTools:::.makeFilter(genoData, c("filt1", "filt2"))
    checkTrue(all(f == c("filt1;filt2", "filt1", "PASS")))
    checkEquals(attributes(f)$header, DataFrame(Description=c("filter 1", "filter 2"), row.names=c("filt1", "filt2")))

    f <- GWASTools:::.makeFilter(genoData, c("filt1"))
    checkTrue(all(f == c("filt1", "filt1", "PASS")))
    checkEquals(attributes(f)$header, DataFrame(Description=c("filter 1"), row.names=c("filt1")))

    f <- GWASTools:::.makeFilter(genoData, NULL)
    checkTrue(all(f == rep("PASS", 3)))
    checkEquals(attributes(f), NULL)
}

test_makeInfo <- function() {
    snp <- data.frame(snpID=1:3,
                      chromosome=1:3,
                      position=1:3,
                      id1=c(TRUE, TRUE, FALSE),
                      id2=c("a", "b", "c"),
                      id3=1:3,
                      id4=c(FALSE,TRUE,TRUE),
                      id5=c(1.0, 1.1, 1.2))
    snpAnnot <- SnpAnnotationDataFrame(snp)
    varMetadata(snpAnnot)[paste0("id", 1:5), "labelDescription"] <-
        paste("id", 1:5)
    mgr <- MatrixGenotypeReader(matrix(0,nrow=3,ncol=2),
                                snpID=snp$snpID,
                                chromosome=snp$chromosome,
                                position=snp$position,
                                scanID=1:2)
    genoData <- GenotypeData(mgr, snpAnnot=snpAnnot)

    i <- GWASTools:::.makeInfo(genoData, paste0("id", 1:5))
    checkTrue(all(all(DataFrame(snp[,4:8]) == i)))
    checkEquals(metadata(i)$header, DataFrame(Number=c(0,1,1,0,1),
                                              Type=c("Flag", "String", "Integer", "Flag", "Float"),
                                              Description=paste("id",1:5),
                                              row.names=paste0("id",1:5)))

    i <- GWASTools:::.makeInfo(genoData, "id1")
    checkTrue(all(all(DataFrame(snp[,4,drop=FALSE]) == i)))
    checkEquals(metadata(i)$header, DataFrame(Number=0, Type="Flag", Description="id 1", row.names="id1"))

    i <- GWASTools:::.makeInfo(genoData, NULL)
    checkEquals(i, DataFrame(matrix(nrow=3, ncol=0)))
}

.testGenoData <- function(nsnp, nsamp) {
    snp <- data.frame(snpID=1:nsnp,
                      chromosome=1:nsnp,
                      position=1:nsnp,
                      alleleA=rep("A", nsnp),
                      alleleB=rep("G", nsnp),
                      stringsAsFactors=FALSE)
    samp <- data.frame(scanID=1:nsamp)
    mgr <- MatrixGenotypeReader(matrix(2, nrow=nsnp, ncol=nsamp),
                                snpID=snp$snpID,
                                chromosome=snp$chromosome,
                                position=snp$position,
                                scanID=samp$scanID)
    GenotypeData(mgr, snpAnnot=SnpAnnotationDataFrame(snp),
                 scanAnnot=ScanAnnotationDataFrame(samp))
}

test_scan.exclude <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(3,5)
    vcf <- genoDataAsVCF(genoData, scan.exclude=c(2,4))
    checkIdentical(geno(vcf)$GT, matrix("0/0", nrow=3, ncol=3, dimnames=list(1:3, c(1,3,5))))
}

test_snp.exclude <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(5,3)
    vcf <- genoDataAsVCF(genoData, snp.exclude=c(2,4))
    checkIdentical(geno(vcf)$GT, matrix("0/0", nrow=3, ncol=3, dimnames=list(c(1,3,5), 1:3)))
}

test_both.exclude <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(4,4)
    vcf <- genoDataAsVCF(genoData, scan.exclude=1, snp.exclude=1)
    checkIdentical(geno(vcf)$GT, matrix("0/0", nrow=3, ncol=3, dimnames=list(2:4, 2:4)))
}

test_scan.order <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(3,5)
    vcf <- genoDataAsVCF(genoData, scan.order=c(4,3,5,1))
    checkIdentical(colnames(vcf), as.character(c(4,3,5,1)))
    vcf <- genoDataAsVCF(genoData, scan.order=5:1, scan.exclude=2)
    checkIdentical(colnames(vcf), as.character(c(5,4,3,1)))
}

test_ref.allele <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(3,2)
    vcf <- genoDataAsVCF(genoData, ref=rep("A", 3))
    checkIdentical(geno(vcf)$GT, matrix("0/0", nrow=3, ncol=2, dimnames=list(1:3, 1:2)))
    checkIdentical(as.character(ref(vcf)), rep("A", 3))
    checkIdentical(as.character(unlist(alt(vcf))), rep("G", 3))

    vcf <- genoDataAsVCF(genoData, ref=rep("B", 3))
    checkIdentical(geno(vcf)$GT, matrix("1/1", nrow=3, ncol=2, dimnames=list(1:3, 1:2)))
    checkIdentical(as.character(ref(vcf)), rep("G", 3))
    checkIdentical(as.character(unlist(alt(vcf))), rep("A", 3))

    vcf <- genoDataAsVCF(genoData, ref=c("A","B","A"))
    checkIdentical(geno(vcf)$GT, matrix(c("0/0", "1/1", "0/0"), nrow=3, ncol=2, dimnames=list(1:3, 1:2)))
    checkIdentical(as.character(ref(vcf)), c("A","G","A"))
    checkIdentical(as.character(unlist(alt(vcf))), c("G","A","G"))
}

test_rowRanges <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(5,3)
    snp <- getSnpAnnotation(genoData)
    snp$id <- c("rs1", NA, "rs3", NA, "rs5")
    genoData <- GenotypeData(genoData@data, snpAnnot=snp, scanAnnot=genoData@scanAnnot)
    vcf <- genoDataAsVCF(genoData, id.col="id")
    checkIdentical(c("rs1", "2:2_A/G", "rs3", "4:4_A/G", "rs5"), rownames(vcf))
}

## convert vcf to gds and return GenotypeData object
.vcf2gds <- function(vcffile, gdsfile) {
    require(SNPRelate)
    snpgdsVCF2GDS(vcffile, gdsfile, verbose=FALSE)
    gds <- GdsGenotypeReader(gdsfile)
    snp <- data.frame(snpID=getSnpID(gds),
                      chromosome=as.integer(getChromosome(gds)),
                      position=getPosition(gds),
                      alleleA=getAlleleA(gds),
                      alleleB=getAlleleB(gds),
                      rsID=getVariable(gds, "snp.rs.id"),
                      qual=getVariable(gds, "snp.annot/qual"),
                      info="A",
                      stringsAsFactors=FALSE)
    filt <- getVariable(gds, "snp.annot/filter")
    for (v in unique(filt)) {
        snp[[v]] <- filt %in% v
    }
    samp <- data.frame(scanID=getScanID(gds),
                       stringsAsFactors=FALSE)
    GenotypeData(gds, snpAnnot=SnpAnnotationDataFrame(snp),
                 scanAnnot=ScanAnnotationDataFrame(samp))
}

test_snprelate <- function() {
    require(VariantAnnotation)
    origfile <- system.file("extdata", "sequence.vcf", package="SNPRelate")
    gdsfile <- tempfile()
    genoData <- .vcf2gds(origfile, gdsfile)
    newfile <- tempfile()
    vcf <- genoDataAsVCF(genoData, id.col="rsID", qual.col="qual", filter.cols=c("PASS", "q10"), info.cols="info")
    writeVcf(vcf, newfile)
    newgds <- tempfile()
    newGenoData <- .vcf2gds(newfile, newgds)
    checkIdentical(getGenotype(newGenoData), getGenotype(genoData))
    checkIdentical(pData(newGenoData@snpAnnot)[,1:5], pData(genoData@snpAnnot)[,1:5])
    checkIdentical(getScanID(newGenoData), getScanID(genoData))
    close(genoData)
    close(newGenoData)
    unlink(c(gdsfile, newfile, newgds))
}

test_VA <- function() {
    require(VariantAnnotation)
    origfile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
    gdsfile <- tempfile()
    genoData <- .vcf2gds(origfile, gdsfile)
    newvcf <- genoDataAsVCF(genoData, id.col="rsID")

    origvcf <- readVcf(origfile, "hg18",
                       param=ScanVcfParam(geno="GT", info=NA))
    origvcf <- origvcf[rownames(newvcf)]
    
    checkIdentical(geno(newvcf)$GT == "0/0", geno(origvcf)$GT == "0|0")
    checkIdentical(geno(newvcf)$GT == "1/1", geno(origvcf)$GT == "1|1")
    checkIdentical(geno(newvcf)$GT %in% "0/1", geno(origvcf)$GT %in% c("0|1", "1|0"))
    checkIdentical(geno(newvcf)$GT == ".", geno(origvcf)$GT == ".")
    checkIdentical(as.character(ref(newvcf)), as.character(ref(origvcf)))
    checkIdentical(as.character(unlist(alt(newvcf))), as.character(unlist(alt(origvcf))))
    close(genoData)
    unlink(c(gdsfile))
}
