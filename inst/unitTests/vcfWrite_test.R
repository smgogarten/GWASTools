test_setFilter <- function() {
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

    checkIdentical(GWASTools:::.setFilter(genoData, c("filt1", "filt2")),
                   list(filter=c("filt1;filt2", "filt1", "PASS"),
                        meta=c(filt1='##FILTER=<ID=filt1,Description="filter 1">',
                               filt2='##FILTER=<ID=filt2,Description="filter 2">')))

    checkIdentical(GWASTools:::.setFilter(genoData, c("filt1")),
                   list(filter=c("filt1", "filt1", "PASS"),
                        meta=c(filt1='##FILTER=<ID=filt1,Description="filter 1">')))

    checkIdentical(GWASTools:::.setFilter(genoData, NULL),
                   list(filter=rep("PASS", 3),
                        meta=character()))
}

test_setInfo <- function() {
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

    checkIdentical(GWASTools:::.setInfo(genoData, paste0("id", 1:5)),
                   list(info=c("id1;id2=a;id3=1;id5=1",
                            "id1;id2=b;id3=2;id4;id5=1.1",
                            "id2=c;id3=3;id4;id5=1.2"),
                        meta=c(id1='##INFO=<ID=id1,Number=0,Type=Flag,Description="id 1">',
                            id2='##INFO=<ID=id2,Number=1,Type=String,Description="id 2">',
                            id3='##INFO=<ID=id3,Number=1,Type=Integer,Description="id 3">',
                            id4='##INFO=<ID=id4,Number=0,Type=Flag,Description="id 4">',
                            id5='##INFO=<ID=id5,Number=1,Type=Float,Description="id 5">')))

    checkIdentical(GWASTools:::.setInfo(genoData, "id1"),
                   list(info=c("id1", "id1", "."),
                        meta=c(id1='##INFO=<ID=id1,Number=0,Type=Flag,Description="id 1">')))

    checkIdentical(GWASTools:::.setInfo(genoData, NULL),
                   list(info=rep(".", 3),
                        meta=character()))
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
    newfile <- tempfile()
    vcfWrite(genoData, newfile, scan.exclude=c(2,4))
    vcf <- readVcf(newfile, "hg18")
    checkIdentical(geno(vcf)$GT, matrix("0/0", nrow=3, ncol=3, dimnames=list(1:3, c(1,3,5))))
    unlink(newfile)
}

test_snp.exclude <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(5,3)
    newfile <- tempfile()
    vcfWrite(genoData, newfile, snp.exclude=c(2,4))
    vcf <- readVcf(newfile, "hg18")
    checkIdentical(geno(vcf)$GT, matrix("0/0", nrow=3, ncol=3, dimnames=list(c(1,3,5), 1:3)))
    unlink(newfile)
}

test_both.exclude <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(4,4)
    newfile <- tempfile()
    vcfWrite(genoData, newfile, scan.exclude=1, snp.exclude=1)
    vcf <- readVcf(newfile, "hg18")
    checkIdentical(geno(vcf)$GT, matrix("0/0", nrow=3, ncol=3, dimnames=list(2:4, 2:4)))
    unlink(newfile)
}

test_snp.exclude.blocks <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(10,3)
    newfile <- tempfile()
    vcfWrite(genoData, newfile, snp.exclude=c(2,4), block.size=4)
    vcf <- readVcf(newfile, "hg18")
    checkIdentical(geno(vcf)$GT, matrix("0/0", nrow=8, ncol=3, dimnames=list(c(1,3,5:10), 1:3)))
    unlink(newfile)
}

test_ref.allele <- function() {
    require(VariantAnnotation)
    genoData <- .testGenoData(3,2)
    newfile <- tempfile()
    vcfWrite(genoData, newfile, ref=rep("A", 3))
    vcf <- readVcf(newfile, "hg18")
    checkIdentical(geno(vcf)$GT, matrix("0/0", nrow=3, ncol=2, dimnames=list(1:3, 1:2)))
    checkIdentical(as.character(ref(vcf)), rep("A", 3))
    checkIdentical(as.character(unlist(alt(vcf))), rep("G", 3))

    vcfWrite(genoData, newfile, ref=rep("B", 3))
    vcf <- readVcf(newfile, "hg18")
    checkIdentical(geno(vcf)$GT, matrix("1/1", nrow=3, ncol=2, dimnames=list(1:3, 1:2)))
    checkIdentical(as.character(ref(vcf)), rep("G", 3))
    checkIdentical(as.character(unlist(alt(vcf))), rep("A", 3))

    vcfWrite(genoData, newfile, ref=c("A","B","A"))
    vcf <- readVcf(newfile, "hg18")
    checkIdentical(geno(vcf)$GT, matrix(c("0/0", "1/1", "0/0"), nrow=3, ncol=2, dimnames=list(1:3, 1:2)))
    checkIdentical(as.character(ref(vcf)), c("A","G","A"))
    checkIdentical(as.character(unlist(alt(vcf))), c("G","A","G"))
    unlink(newfile)
}

## convert vcf to gds and return GenotypeData object
.vcf2gds <- function(vcffile, gdsfile) {
    require(SNPRelate)
    snpgdsVCF2GDS(vcffile, gdsfile)
    gds <- GdsGenotypeReader(gdsfile)
    snp <- data.frame(snpID=getSnpID(gds),
                      chromosome=as.integer(getChromosome(gds)),
                      position=getPosition(gds),
                      alleleA=getAlleleA(gds),
                      alleleB=getAlleleB(gds),
                      rsID=getVariable(gds, "snp.rs.id"),
                      stringsAsFactors=FALSE)
    samp <- data.frame(scanID=getScanID(gds),
                       stringsAsFactors=FALSE)
    GenotypeData(gds, snpAnnot=SnpAnnotationDataFrame(snp),
                 scanAnnot=ScanAnnotationDataFrame(samp))
}

test_snprelate <- function() {
    origfile <- system.file("extdata", "sequence.vcf", package="SNPRelate")
    gdsfile <- tempfile()
    genoData <- .vcf2gds(origfile, gdsfile)
    newfile <- tempfile()
    vcfWrite(genoData, newfile, id.col="rsID")
    newgds <- tempfile()
    newGenoData <- .vcf2gds(newfile, newgds)
    checkIdentical(getGenotype(newGenoData), getGenotype(genoData))
    checkIdentical(pData(newGenoData@snpAnnot), pData(genoData@snpAnnot))
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
    newfile <- tempfile()
    vcfWrite(genoData, newfile, id.col="rsID")
    newvcf <- readVcf(newfile, "hg18")

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
    unlink(c(gdsfile, newfile))
}


test_vcfCheck <- function() {
    genoData <- .testGenoData(5,3)
    newfile <- tempfile()
    vcfWrite(genoData, newfile)
    vcfCheck(genoData, newfile)
    unlink(newfile)
}

test_vcfCheck_ref <- function() {
    genoData <- .testGenoData(5,3)
    newfile <- tempfile()
    vcfWrite(genoData, newfile, ref.allele=c("A","B","A","B","A"))
    vcfCheck(genoData, newfile)
    unlink(newfile)
}
