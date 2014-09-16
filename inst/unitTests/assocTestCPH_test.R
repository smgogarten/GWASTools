## fixing a bug for data with no Y chromosome
test_setSCF <- function() {
    scf <- GWASTools:::.setSCF(NULL, scanID=1:5, chrom=rep(1:24, each=2), ychrom.code=25,
                   sex=rep("F", 5), cvnames="sex", scan.exclude=NULL)
    checkIdentical(matrix(TRUE, nrow=5, ncol=24, dimnames=list(1:5, 1:24)), scf)

    scf <- GWASTools:::.setSCF(NULL, scanID=1:5, chrom=rep(1:25, each=2), ychrom.code=25,
                   sex=rep("F", 5), cvnames="sex", scan.exclude=NULL)
    checkIdentical(matrix(c(rep(TRUE, 5*24), rep(FALSE, 5)),
                          nrow=5, ncol=25, dimnames=list(1:5, 1:25)), scf)

    scf <- GWASTools:::.setSCF(NULL, scanID=1:5, chrom=rep(1:25, each=2), ychrom.code=25,
                   sex=rep("M", 5), cvnames="sex", scan.exclude=NULL)
    checkIdentical(matrix(c(rep(TRUE, 5*24), rep(FALSE, 5)),
                          nrow=5, ncol=25, dimnames=list(1:5, 1:25)), scf)

    scf <- GWASTools:::.setSCF(NULL, scanID=1:5, chrom=rep(1:24, each=2), ychrom.code=25,
                   sex=rep("F", 5), cvnames="sex", scan.exclude=4:5)
    checkIdentical(matrix(c(TRUE, TRUE, TRUE, FALSE, FALSE),
                          nrow=5, ncol=24, dimnames=list(1:5, 1:24)), scf)
}
