#### Circular Binary Segmentation #############


########### Main function ####################
anomSegmentBAF<-function(intenData, genoData, scan.ids, chrom.ids, snp.ids,
                         smooth=50, min.width=5, nperm=10000, alpha=.001,
                         verbose=TRUE){
#scan.ids: vector of sample numbers to process
#chrom.ids: vector of (unique) chromosomes to process
#snp.ids: vector of eligible snp ids
## smooth: number of markers for smoothing region for DNAcopy
## min.width: minimum number of markers for segmenting in DNAcopy
## nperm: number of permutations for DNAcopy
## alpha: significance level for DNAcopy
## verbose: parameter to DNAcopy 0=no output, 1=prints sample number currently processing

  # check that intenData has BAF
  if (!hasBAlleleFreq(intenData)) stop("BAlleleFreq not found in intenData")
  
  # check that dimensions of intenData and genoData are equal
  intenSnpID <- getSnpID(intenData)
  genoSnpID <- getSnpID(genoData)
  if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
  intenScanID <- getScanID(intenData)
  genoScanID <- getScanID(genoData)
  if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")
  
  intid <- intenSnpID
  if(!all(is.element(snp.ids,intid))) stop("eligible snps not contained in snp ids")

  anoms<-NULL
  ## compute parameter needed for DNAcopy
  max.ones<-floor(nperm*alpha)+1
  sbdry<-getbdry(.05,nperm,max.ones)

  sid <- intenScanID
  chrom <- getChromosome(intenData)
  for(snum in scan.ids){
    sindex <- which(is.element(sid, snum))
    if(length(sindex)==0) stop(paste("Sample ",snum, " does not exist",sep=""))
    GENO <- getGenotype(genoData, snp=c(1,-1), scan=c(sindex,1))
  # get genotypes for the given sample for all snps

  ## consider only hets (1) and missing: intensity only and failed snps excluded
    sel<-is.element(intid,snp.ids) & (GENO == 1 | is.na(GENO))

    baf <- getBAlleleFreq(intenData, snp=c(1,-1), scan=c(sindex,1))
    ws<-!is.na(baf)

    INDEX<-which(sel&ws)
    CHR<-chrom[sel&ws]
    BAF<-baf[sel&ws]

    index<-NULL
    chr<-NULL
    baf.dat<-NULL
    for(ch in chrom.ids){
      wc<-CHR==ch    #T/F for indices of CHR which match indices of INDEX,BAF
      bf<-BAF[wc]
      if(length(bf)<3) next  # don't want to feed DNAcopy something that will crash
      ind<-INDEX[wc]
      chrm<-CHR[wc]
      med<-median(bf,na.rm=T)
      bf1<-1-bf
      bfm<-abs(bf-med)
      c<-cbind(bf,bf1,bfm)
      met<-apply(c,1,min)
      baf.metric<-sqrt(met)
      index<-c(index,ind)
      chr<-c(chr,chrm)
      baf.dat<-c(baf.dat,baf.metric)   
    } #end of chrom loop

    temp.CNA<-CNA(as.vector(baf.dat),chr,index,data.type="logratio",sampleid=snum)
    temp.smooth<-smooth.CNA(temp.CNA,smooth.region=smooth,outlier.SD.scale=4)  #smooth.region,outlier.SD.scale=default
    temp.segment<-segment(temp.smooth,alpha=alpha,sbdry=sbdry,p.method="h",min.width=min.width,nperm=nperm,undo.splits="sdundo",undo.SD=1,verbose=as.integer(verbose))
    tmp<-temp.segment$out
    if(dim(tmp)[1]<1) next
    tmp$ID<-snum
    tmp$loc.start<-as.integer(tmp$loc.start)
    tmp$loc.end<-as.integer(tmp$loc.end)
    #$loc.start and $loc.end are indices of snp/indices of intid

    names(tmp)<-c("scanID","chromosome","left.index","right.index","num.mark","seg.mean")
    anoms<-rbind(anoms,tmp)
  } #end loop on samples

  
  return(anoms)
}#end function
