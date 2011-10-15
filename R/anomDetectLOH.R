
anomDetectLOH<-function(intenData, genoData, scan.ids, chrom.ids, snp.ids,
                        known.anoms,
                        smooth=50,min.width=5,nperm=10000,alpha=.001,
                        run.size=50,inter.size=4,
                        homodel.min.num=10,homodel.thresh=10,
                        small.num=20,small.thresh=2.25, medium.num=50,medium.thresh=2,
                        long.num=100,long.thresh=1.5, small.na.thresh=2.5,
                        length.factor=5,merge.fac=.85,min.lrr.num=20,verbose=TRUE){

##scan.ids: samples to consider
##chrom.ids: chromosomes to consider, usually 1:23
##snp.ids: intid's for eligible snps
##known.anoms: data.frame of known anoms (usually from BAF); 
  #   must have "scanID","chrom","left","right" where left and right are snp indices

##run.size: number of markers to declare a 'homozygous' run (here homozygous includes missing) 
##inter.size: number of consecutive heterozygous markers allowed to interrupt a run
### detection of anomalies based on a chromsome-wide and local mad.fac thresholds
#    where mad.fac is (segment median-nonanomalous median)/nonanom mad
##homodel.min.num: minimum number of markers to detect extreme diff in lrr (for homo deletion)
##homodel.thresh: threshold for 'mad.fac' to detect extreme diff in lrr
##small.num: minimum number of markers to detect anomaly (other than extreme)
##small.thresh:threshold for  number of markers between small.num and medium.num
##medium.num: number of markers threshold for 'medium' 
##medium.thresh: threshold for 'mad.fac' when number of markers is between medium num and long num
##long.num: number of markers threshold for 'long'
##long.thresh:threshold for 'mad.fac' when number of markers is bigger than long.num
##small.na.thresh: 
#   chrom mad.fac threshold when between small.num and medium.num and no local mad.fac
##length.factor: 
#   local mad.fac based on interval that is length.factor*(no. of markers in segment) on either side
##merge.fac: threshold used to merge original segmentation segments 
##min.lrr.num: if any 'nonanomalous' interval less than min.lrr.num,
#  ignore this piece in finding overall nonanomalous unless is only piece left

#### checks ####
  # check that intenData has LRR
  if (!hasLogRRatio(intenData)) stop("LogRRatio not found in intenData")
  
  # check that dimensions of intenData and genoData are equal
  intenSnpID <- getSnpID(intenData)
  genoSnpID <- getSnpID(genoData)
  if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
  intenScanID <- getScanID(intenData)
  genoScanID <- getScanID(genoData)
  if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")
  
  # check that sex is present in annotation
  if (hasSex(intenData)) {
    sex <- getSex(intenData)
  } else if (hasSex(genoData)) {
    sex <- getSex(genoData)
  } else stop("sex not found in intenData or genoData")

  intid <- intenSnpID
  if(!all(is.element(snp.ids,intid))) stop("eligible snps not contained in snp ids")

  sid <- intenScanID
  male <- sid[sex == "M"]
  chrom <- getChromosome(intenData)
  
  if(!is.element(class(known.anoms),"data.frame") | !all(is.element(c("scanID","chromosome","left.index","right.index"),names(known.anoms)))){
    stop("known.anoms input needs to be data.frame with variables including scanID, chromosome, left.index, right.index") }
####

  # internal functions require these names, so convert from package standard
  names(known.anoms)[names(known.anoms) == "chromosome"] <- "chrom"
  names(known.anoms)[names(known.anoms) == "left.index"] <- "left"
  names(known.anoms)[names(known.anoms) == "right.index"] <- "right"

  LOH.raw<-NULL;LOH.base.info<-NULL; LOH.filtered<-NULL
  LOH.segments<-NULL; LOH.merge<-NULL;LOH.raw.adjusted<-NULL
  ## compute parameter needed for DNAcopy
  max.ones<-floor(nperm*alpha)+1
  sbdry<-getbdry(.05,nperm,max.ones)

  sel<-is.element(intid,snp.ids)
  orindex<-which(sel)  #indices of eligible snp's
  orchr<-chrom[sel]
  for(snum in scan.ids){
    sindex <- which(is.element(sid, snum))
    if(length(sindex)==0) stop(paste("Sample ",snum, " does not exist",sep=""))
    GENO <- getGenotype(genoData, snp=c(1,-1), scan=c(sindex,1))
    olrr<-getLogRRatio(intenData, snp=c(1,-1), scan=c(sindex,1))[sel]  
    ogeno<-GENO[sel]
    ws<-!is.na(olrr)
    olrr<-olrr[ws]
    oindex<-orindex[ws]
    ogeno<-ogeno[ws]
    chr<-orchr[ws]

    ##homoz and missing have value 0
    whomo<-is.element(ogeno,c(0,2)) | is.na(ogeno)
    ogeno[whomo]<-0

    ####### DNAcopy Segmentation for given sample ########
    temp.CNA<-CNA(as.vector(olrr),chr,oindex,data.type="logratio",sampleid=snum)
    temp.smooth<-smooth.CNA(temp.CNA,smooth.region=smooth,outlier.SD.scale=4)  #smooth.region,outlier.SD.scale=default
    temp.segment<-segment(temp.smooth,alpha=alpha,sbdry=sbdry,p.method="h",min.width=min.width,nperm=nperm,undo.splits="sdundo",undo.SD=1,verbose=as.integer(verbose))
    segments<-temp.segment$out
    if(dim(segments)[1]<1) next
    segments$ID<-as.integer(rep(snum,length(segments$ID)))
    #$loc.start and $loc.end are indices of snp
    names(segments)<-c("scanID","chrom","left","right","num.mark","seg.mean")


    for(ch in chrom.ids){
      if(ch==XchromCode(intenData) & is.element(snum,male)) next

      wc<-chr==ch 
      geno<-ogeno[wc]
      index<-oindex[wc]
      lrr<-olrr[wc]   
      chrr<-chr[wc]

      ansch<-known.anoms[known.anoms$scanID==snum & known.anoms$chrom==ch,]
      segs<-segments[segments$chrom==ch,]

      ###### find homozygous runs and base info for current sample/chrom
      out<-GWASTools:::LOHfind(snum,ch,geno,index,lrr,chrr,segs,ansch,
                   run.size,inter.size,min.lrr.num )
  
      base.snch<-out$base.info
      base.snch$sex<-sex[sindex]
      RUNS.snch<-out$RUNS
      segs.snch<-out$segments

      LOH.base.info<-rbind(LOH.base.info,base.snch)
      LOH.segments<-rbind(LOH.segments,segs.snch)
      LOH.raw<-rbind(LOH.raw,RUNS.snch)

      #### BEGIN filtering process  #############
  
       ####  Special Cases ###########  
    
      if(base.snch$num.runs==0) next  #no anomalies

      if(is.na(base.snch$chrom.nonanom.mad) & is.element(base.snch$num.runs,c(1,2)) & dim(ansch)[1]==0) {
        # happens if runs take up most of chromosome with no other known anoms, i.e. probably whole chrom LOH
        seg.median<-NULL;seg.mean<-NULL;nm<-NULL
        for(k in 1:dim(RUNS.snch)[1]){
          subind<-index>=RUNS.snch$left[k] & index<=RUNS.snch$right[k]         
          nm<-c(nm,sum(subind))
          sublrr<-lrr[subind]
          seg.med<-median(sublrr,na.rm=T)
          seg.median<-c(seg.median,seg.med)
          seg.mean<-c(seg.mean,mean(sublrr,na.rm=T))
        }#end of k loop
        RUNS.snch$num.mark<-nm
        RUNS.snch$seg.median<-seg.median
        RUNS.snch$seg.mean<-seg.mean
        RUNS.snch$sd.fac<-NA
        RUNS.snch$mad.fac<-NA
        RUNS.snch$local<-NA
        RUNS.snch$num.segs<-base.snch$num.segs
        RUNS.snch$chrom.nonanom.mad<-NA
        RUNS.snch$chrom.nonanom.median<-NA
        RUNS.snch$chrom.nonanom.mean<-NA
        RUNS.snch$chrom.nonanom.sd<-NA
        RUNS.snch$sex<-sex[sindex] 
        LOH.raw.adjusted<-rbind(LOH.raw.adjusted,RUNS.snch)
        LOH.filtered<-rbind(LOH.filtered,RUNS.snch)
        next
      }
  
      if(base.snch$num.segs==1) next
      #if no segmentation and previous situation not occur, there are no anomalies
      if(is.na(base.snch$chrom.nonanom.mad)) next  
      #no base left to compare with;already singled out whole chrom possibiltiy; too many pieces: no anoms

      ##### end Special cases###
      ### Get info needed for filtering process
      ## find indices for 'nonanom' snps: not in known anoms nor in any found runs (which are potential anoms) 
      baf.del<-NULL
      if(dim(ansch)[1]!=0){
        for(j in 1:dim(ansch)[1]) {
          int<-index[index>=ansch$left[j]& index<=ansch$right[j]]
          baf.del<-union(baf.del,int)
        }
      }  #index values
  
      runs.del<-NULL             #delete all found runs or anoms to determine "nonanom" base
      for(i in 1:dim(RUNS.snch)[1]){
        int<-index[index>=RUNS.snch$left[i] & index<=RUNS.snch$right[i]]
        runs.del<-union(runs.del,int) 
      }    
  
      possible.anom.index<-union(baf.del,runs.del)
      #by previous checks, will have at least some runs and runs+known anoms will not take up whole chrom
      nonanom.index<-setdiff(index,possible.anom.index)

      #### final FILTERING #########
      outt<-GWASTools:::LOHselectAnoms(snum,ch,segs.snch,RUNS.snch,base.snch,index,nonanom.index,lrr,
        homodel.min.num,homodel.thresh,small.num,small.thresh, medium.num,medium.thresh,
        long.num,long.thresh, small.na.thresh,length.factor,merge.fac,min.lrr.num)
      raw.adj<-outt$raw.adjusted
      if(!is.null(raw.adj)) raw.adj$sex<-sex[sindex]  
      filtered<-outt$filtered
      if(!is.null(filtered)) filtered$sex<-sex[sindex] 
      LOH.raw.adjusted<-rbind(LOH.raw.adjusted,raw.adj)
      LOH.filtered<-rbind(LOH.filtered, filtered)
      if(outt$merge.flag){
        tmp<-data.frame(snum,ch);names(tmp)<-c("scanID","chrom")
        LOH.merge<-rbind(LOH.merge,tmp)
      }

    } #end chrom loop
  } #end sample loop


  colm<-c("scanID","chrom","left","right","num.mark","seg.median","seg.mean",
          "mad.fac","sd.fac","local","num.segs","chrom.nonanom.mad","chrom.nonanom.median","chrom.nonanom.mean","chrom.nonanom.sd","sex")
  
  if(!is.null(LOH.raw.adjusted)){
    LOH.raw.adjusted<-LOH.raw.adjusted[,colm]
  }

  if(!is.null(LOH.filtered)){
    LOH.filtered<-LOH.filtered[order(LOH.filtered$scanID,LOH.filtered$chrom,LOH.filtered$left),]
    LOH.filtered<-LOH.filtered[,colm] 
  }

  outr<-list(LOH.raw,LOH.raw.adjusted,LOH.filtered,LOH.base.info,LOH.segments,LOH.merge)
  names(outr)<-c("raw","raw.adjusted","filtered","base.info","segments","merge")

  # convert back to package standard names
  for (i in 1:length(outr)) {
    if ("chrom" %in% names(outr[[i]]))
      names(outr[[i]])[names(outr[[i]]) == "chrom"] <- "chromosome"
    if ("left" %in% names(outr[[i]]))
      names(outr[[i]])[names(outr[[i]]) == "left"] <- "left.index"
    if ("right" %in% names(outr[[i]]))
      names(outr[[i]])[names(outr[[i]]) == "right"] <- "right.index"
  }

  # convert index to base position
  for (i in 1:length(outr)) {
    if (!is.null(outr[[i]]) & ("left.index" %in% names(outr[[i]]))) {
      outr[[i]]$left.base <- getPosition(intenData, index=outr[[i]]$left.index)
      outr[[i]]$right.base <- getPosition(intenData, index=outr[[i]]$right.index)
    }
  }
    
  return(outr)
}


