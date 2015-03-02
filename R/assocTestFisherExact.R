assocTestFisherExact <- function(dat, outfile = NULL){

    .Deprecated("batchFisherTest")
    
  # columns that contain count data
  count.cols <- grep(paste("n[[:upper:]]{2}", sep="."), names(dat))
  # extract genotype counts
  # 6 columns:  nAA.cc0, nAB.cc0, nBB.cc0, nAA.cc1, nAB.cc1, nBB.cc1
  count.dat <- dat[,count.cols]

  # number of SNPs
  nsnp <- dim(count.dat)[1]

  # change genotype counts to allele counts
  allele.dat <- matrix(NA, nrow=nsnp, ncol=4)
  allele.dat[,1] <- 2*count.dat[,1]+count.dat[,2]
  allele.dat[,2] <- 2*count.dat[,3]+count.dat[,2]
  allele.dat[,3] <- 2*count.dat[,4]+count.dat[,5]
  allele.dat[,4] <- 2*count.dat[,6]+count.dat[,5]
  allele.dat <- as.data.frame(allele.dat)
  names(allele.dat) <- c("nA.cc0", "nB.cc0", "nA.cc1", "nB.cc1")


  # index for column containing minor allele info
  maidx <- grep("minor.allele", names(dat))
  
  # results matrix
  res <- matrix(NA, nrow=nsnp, ncol=4)
    
  # perform Fisher's Exact test
  for(i in 1:nsnp){
    # if all allele counts 0
    if(is.na(dat[i,maidx])){
      res[i,] <- rep(NA, 4)
    }else{
      # if minor allele is A
      if(dat[i,maidx]=="A"){
        tmp <- matrix(as.numeric(allele.dat[i,c(3,4,1,2)]), nrow=2)
      # if minor allele is B
      }else if(dat[i,maidx]=="B"){
        tmp <- matrix(as.numeric(allele.dat[i,c(4,3,2,1)]), nrow=2)
      }
      res[i,] <- as.numeric(unlist(fisher.test(tmp))[c(4,2,3,1)])
    }
  }
  res <- as.data.frame(res)
  names(res) <- c("Fisher.OR","Fisher.OR_L95","Fisher.OR_U95","Fisher.pval")

  # index for column containing MAF info
  mafidx <- grep("MAF", names(dat))
  # index for column containing warningOrError info
  weidx <- grep("warningOrError", names(dat))
 
  # combine results
  res <- cbind(dat[c(1,2,mafidx,maidx,weidx)], res, allele.dat)
  names(res)[c(2:5)] <- c("n","MAF","minor.allele","regression.warningOrError")

  # save results
  if(!is.null(outfile)){
    fileOut <- paste(outfile, "FisherExact", "RData", sep=".")
    save(res, file=fileOut, compress=TRUE)
  }else{
    return(res)
  }
}
