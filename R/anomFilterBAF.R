
#############SUB-FUNCTIONS############
## MERGE contiguous segments that meet a filter
#  based on number of st.devs from a baseline
mergeSeg<-function(segs,snum,ch,cent,Pos,base.mean,base.sd,sd.reg,sd.long,num.mark.thresh,long.num.mark.thresh,low.frac.used) { 
#segs is data.frame of DNAcopy segments from sample snum and chromosome ch
#cent: centromere position information
#Pos: position info from netCDF
if(!is.element(class(segs),"data.frame")) stop("segment data is not a data.frame")
if(dim(segs)[1]==0)stop("error: no segment info given to merge")
colm<-c("scanID","chrom","left","right","num.mark","seg.mean","sd.fac","sex")
if(!all(is.element(colm,names(segs)))){
 stop("incorrect column variables for segs - must be c(scanID,chrom,left,right,num.mark,seg.mean,sd.fac,sex)")}

if(!all(segs$scanID==snum))stop("data.frame of segments needs to be from same sample")
if(!all(segs$chrom==ch))stop("data.frame of segments needs to be from same chromosome")

sx<-segs$sex[1]
an<-segs[order(segs$left),]
an$merge<-FALSE
if(dim(an)[1]<2) return(an)
frac.used<-an$num.mark/(an$right-an$left+1)
  # denom is total number of markers in between, including intensity only
  # this would make it more likely for frac.used to be smaller
  # would be an unusual anom and would probably not want to merge
  
  merged.anoms<-NULL

  s1<-an$sd.fac>=sd.reg
  s2<-an$num.mark>long.num.mark.thresh & an$sd.fac>=sd.long
  s<-s1|s2
  s3<-frac.used>low.frac.used
  ss<-s & s3  # T for ones above threshold and not low.frac
  r<-rle(ss)
  vals<-r[[2]]
  lens<-r[[1]]
  nv<-length(vals)
  endp<-cumsum(lens)  #end positions of each run
  stp<-c(1,endp[1:(nv-1)]+1)
  wt<-which(vals & lens>=2)
  if(length(wt)==0) return(an)
  
  endt<-endp[wt]
  stt<-stp[wt]
  
  del.merge<-NULL
  for(i in 1:length(stt)){
    ind<-stt[i]:endt[i]
    tmp<-an[ind,]  #set of consecutive anoms
    
    ## don't want to merge any anoms with num.mark<5 on either 'edge'
    ## and do not merge if would create a centromere spanning anom
    w5<-which(tmp$num.mark>=5)
    if(length(w5)<=1) next
    mw5<-min(w5);mx5<-max(w5) 
    choose<-mw5:mx5
    nt<-length(choose)
    if(nt<2) next
    if(nt==2) { 
       c<-Pos[tmp[mw5,"right"]]<=cent$left.base[cent$chrom==ch] &
          Pos[tmp[mx5,"left"]]>=cent$right.base[cent$chrom==ch] 
       nm<-tmp[mw5,"num.mark"]>=num.mark.thresh & tmp[mx5,"num.mark"]>=num.mark.thresh
         if(!c & !nm) {
           del.merge<-c(del.merge,ind[choose])
           new.left<-an[ind[mw5],"left"]; new.right<-an[ind[mx5],"right"]
           new.num.mark<-sum(an$num.mark[ind[choose]])
           new.seg.mean<-sum(an$seg.mean[ind[choose]]*an$num.mark[ind[choose]])/new.num.mark
           new.sdfac<-abs(new.seg.mean-base.mean)/base.sd
           new<-data.frame(snum,ch,new.left,new.right,new.num.mark,new.seg.mean,new.sdfac,sx,T,stringsAsFactors=FALSE)
           names(new)<-c("scanID","chrom","left","right","num.mark","seg.mean","sd.fac","sex","merge")
           merged.anoms<-rbind(merged.anoms,new)
         }
      } else {
        I<-mx5
        for(j in mw5:(mx5-1)){
          if(Pos[an[ind[j],"right"]]<=cent$left.base[cent$chrom==ch] &
          Pos[an[ind[j+1],"left"]]>=cent$right.base[cent$chrom==ch]) {I<-j;break} 
        }
        
        set1<-mw5:I ; if(I<mx5) set2<-(I+1):mx5 else set2<-NULL
        if(length(set1)>=2) { del.merge<-c(del.merge,ind[set1])
           new.left<-an[ind[mw5],"left"]; new.right<-an[ind[I],"right"]
           new.num.mark<-sum(an$num.mark[ind[set1]])
           new.seg.mean<-sum(an$seg.mean[ind[set1]]*an$num.mark[ind[set1]])/new.num.mark
           new.sdfac<-abs(new.seg.mean-base.mean)/base.sd
           new<-data.frame(snum,ch,new.left,new.right,new.num.mark,new.seg.mean,new.sdfac,sx,T,stringsAsFactors=FALSE)
           names(new)<-c("scanID","chrom","left","right","num.mark","seg.mean","sd.fac","sex","merge")
           merged.anoms<-rbind(merged.anoms,new)
         }
         if(length(set2)>=2) { del.merge<-c(del.merge,ind[set2])
           new.left<-an[ind[I+1],"left"]; new.right<-an[ind[mx5],"right"]
           new.num.mark<-sum(an$num.mark[ind[set2]])
           new.seg.mean<-sum(an$seg.mean[ind[set2]]*an$num.mark[ind[set2]])/new.num.mark
           new.sdfac<-abs(new.seg.mean-base.mean)/base.sd
           new<-data.frame(snum,ch,new.left,new.right,new.num.mark,new.seg.mean,new.sdfac,sx,T,stringsAsFactors=FALSE)
           names(new)<-c("scanID","chrom","left","right","num.mark","seg.mean","sd.fac","sex","merge")
           merged.anoms<-rbind(merged.anoms,new)
         }
      }#end of else
 } #end of i loop

   if(length(del.merge)!=0){
      tmp<-an[-del.merge,]
      tmpn<-merged.anoms
      out<-rbind(tmp,tmpn)
      out<-out[order(out$left),]
   } else out<-an
return(out)
} #end function
################################

###### delHomoRuns ##############
## function to possibly narrow segments found containing homo del
# look for adjustment for selected anoms
# to identify homozygous deletions
# looking for run of lrr values < lrr.cut then narrow to this run
# (BAF DNAcopy tends to not segment these well - often occur in longer homozygous runs)

delHomoRuns<-function(anoms,sid,eligible,intid,LRR,run.size,inter.size,
   low.frac.used,lrr.cut,ct.thresh,frac.thresh){
#run.size - min length of run
#inter.size - number of homozygotes allowed to "interrupt" run
#low.frac.used - fraction of markers used compared to number of markers in interval
#lrr.cut-look for runs of lrr values below lrr.cut
#ct.thresh - minimum number of lrr values below lrr.cut needed in order to process
#frac.thresh - process only if (# lrr values below lrr.cut)/(# eligible lrr in interval) > frac.thresh
#anoms: data.frame of anomalies  scanID, chrom, left, right, num.mark, seg.mean, sd.fac, sex, merge

if(!is.element("data.frame",class(anoms))) stop("anoms needs to be a data.frame")
annames<-c("scanID","chrom","left","right","num.mark","seg.mean","sd.fac","sex","merge")
if(!all(is.element(annames,names(anoms)))) stop("anoms does not have required variable names")
if(length(unique(anoms$scanID))!=1) stop("anoms needs to be for one sample")

anoms.rev<-NULL

for(I in 1:dim(anoms)[1]){
  an<-anoms[I,]
  snum<-an$scanID; chr<-an$chrom
  ledge<-an$left;redge<-an$right
  sindex <- which(is.element(sid, snum))
  if(length(sindex)==0) stop(paste("Sample ",snum, " does not exist",sep=""))

##want to look for runs only in the already identified anomaly
  int<-intid>=intid[ledge] & intid<=intid[redge]  #$left and $right are indices of intid
  selgood<-is.element(intid,eligible)
  index<-which(selgood&int)

  lrr <- LRR[index]
  wn<-!is.na(lrr)
  lrr<-lrr[wn]
  index<-index[wn]

  whm<-lrr< lrr.cut
  ct<-sum(whm)
  frac<-ct/length(index)
  an$low.ct<-ct
  an$nmark.lrr<-length(index)
  pct<-an$num.mark/(an$right-an$left+1)
  w<-frac>frac.thresh | pct<low.frac.used
  if(ct< ct.thresh |!w){
     an$old.left<-an$left; an$old.right<-an$right
     anoms.rev<-rbind(anoms.rev,an)
     next
  }
  who<-!whm
  lrr[whm]<-0
  lrr[who]<-1

  rgt<-NULL;lft<-NULL
  w<-rle(as.vector(lrr))
  #w<-rle(  ) has w[[1]] = lengths of runs, w[[2]] corresponding values
  vals<-w[[2]];lngs<-w[[1]]
  r0<-vals

  rlen<-length(r0)
  if(rlen==1){  # keep original if only one value    
     an$old.left<-an$left; an$old.right<-an$right
     anoms.rev<-rbind(anoms.rev,an)
     next  
  }

##establish initial positions of alternating runs
  endp<-cumsum(lngs)
  in.pos<-c(1,endp[1:(rlen-1)]+1)

##merging intervals if separated by < inter.size no. of undesirable values
# assuming sum of lengths of desirable intervals on either side meets run.size criterion

  tpos<-which(r0==1&lngs<=inter.size) #identify small runs of homos
  smf<-which(r0==0 & lngs<run.size) #identify 'small' runs of hets
  if(length(tpos)!=0){
    if(tpos[1]==1) {
     if(lngs[2]>=run.size){r0[1]<-0}
     tpos<-tpos[-1]}
     if(length(tpos)!=0){
       if(tpos[length(tpos)]==rlen) {
          tpos<-tpos[-length(tpos)]; if(lngs[rlen-1]>=run.size) {r0[rlen]<-0}}
     if(length(tpos)!=0){
        for(k in tpos){if((lngs[k-1]+lngs[k+1])>=run.size) {r0[k]<-0}
        }
  }} }

##want smaller runs of 0s to become runs of undesirable but not if they
 #are part of a combined run
run2<-rle(r0)
vals2<-run2[[2]];lngs2<-run2[[1]]
w0<-which(vals2==0 & lngs2>1)
if(length(w0)==0){r0[smf]<-1} else {
index.set<-NULL
for(j in 1:length(w0)){ if(w0[j]==1){start<-1} else {start<-sum(lngs2[1:(w0[j]-1)])+1}
end<-sum(lngs2[1:w0[j]])
index.set<-c(index.set,start:end)}
smf.use<-setdiff(smf,intersect(smf,index.set))
r0[smf.use]<-1}

##after merging some of the initial runs, we get modified listing of r0 vals
## look for runs here; e.g. if run of two 0s that means that initially
#we had a run of desirable with small run of undesirable - rle length of 2 then indicates
#putting those two original runs together as one run of desirable 
new.rle<-rle(r0)
nvals<-new.rle[[2]]
nlens<-new.rle[[1]]  ##indicates how many of original runs to put together
if(length(nvals)==1){ #all now classified as undesirable or as desirable = no change
   an$old.left<-an$left; an$old.right<-an$right
   anoms.rev<-rbind(anoms.rev,an)
   next  
 }
 
newt<-which(nvals==0) #newt could be empty if originally there were no long het/miss runs
if(length(newt)==0){     
    an$old.left<-an$left; an$old.right<-an$right
    anoms.rev<-rbind(anoms.rev,an)
    next  
 }
 

left<-NULL
right<-NULL
### if newt indicates runs of 0s, change initial and end positions of 
#runs of desired accordingly
if(newt[1]==1){left<-c(left,1);right<-c(right,in.pos[nlens[1]+1]-1);newt<-newt[-1]}
for(k in newt){ind<-sum(nlens[1:(k-1)]);left<-c(left,in.pos[ind+1])
kl<-length(newt);kk<-newt[kl]
if((ind+1+nlens[kk])<=length(in.pos)){
right<-c(right,in.pos[ind+1+nlens[k]]-1)} else {right<-c(right,length(index))}}
##right and left positions are indices of lrr (= indices of index)

## if splits into more than one run, leave as the original 
if(length(right)==0|length(left)==0|length(right)>1|length(left)>1){
     an$old.left<-an$left; an$old.right<-an$right
     anoms.rev<-rbind(anoms.rev,an)
     next  
   }

## there is one adjusted interval found

an$old.left<-an$left;an$old.right<-an$right
an$left<-index[left];an$right<-index[right] 
anoms.rev<-rbind(anoms.rev,an)
 } #end loop on anomalies
return(anoms.rev)
 }
###################################################
########## MAIN FUNCTION ##########################
##########
anomFilterBAF<-function(intenData, genoData, segments, snp.ids,
   centromere,low.qual.ids=NULL,
  num.mark.thresh=15,long.num.mark.thresh=200,sd.reg=2,sd.long=1,
  low.frac.used=.1,run.size=10,inter.size=2,low.frac.used.num.mark=30,
   very.low.frac.used=.01,low.qual.frac.num.mark=150,
  lrr.cut= -2,ct.thresh=10,frac.thresh=.1,verbose=TRUE){
##segments - data.frame to determine which are anomalous
 # names of data.frame must include "scanID","chromosome","num.mark","left.index","right.index","seg.mean"
 # assume have segmented each chromosome (at least all autosomes) for a given sample
##snp.ids: vector of eligible snp ids
##centromere: data.frame with centromere position info
## low.qual.ids: sample numbers determined to be messy for which segments are filtered
#     based on num.mark and fraction used
## num.mark.thresh: minimum size of segment to consider for anomaly
## long.num.mark.thresh: min number of markers for "long" segment
 # (significance threshold allowed to be lower)
## sd.reg: number of standard deviations from segment mean 
 # compared to a baseline mean for "normal" needed to declare segment anomalous
## sd.long: same meaning as sd. long but applied to "long" segments
## low.frac.used: fraction used to declare a segment with 
 # low number hets or missing compared with number of markers in interval 
## run.size, inter.size: for possible determination of homozygous deletions 
 # (see description in delHomoRuns function above)
## low.frac.used.num.mark: used in final step of deleting "small" low frac.used segments
 # which tend to be false positives (after determining homo deletions)
## very.low.frac.used: any segments with (num.mark)/(number of markers in interval) less than this are filtered out
## low.qual.frac.num.mark: num.mark threshold for messy samples
##lrr.cut-look for runs of lrr values below lrr.cut to adjust homo del endpts
##ct.thresh - minimum number of lrr values below lrr.cut needed in order to adjust
##frac.thresh - adjust only if (# lrr values below lrr.cut)/(# eligible lrr in interval) > frac.thresh

  # check that intenData has BAF
  if (!hasBAlleleFreq(intenData)) stop("BAlleleFreq not found in intenData")
  
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
  chrom <- getChromosome(intenData)
  Pos <- getPosition(intenData)
  
  if(!is.element(class(segments),"data.frame")) stop("data is not a data.frame")
  chk<-is.element(c("scanID","chromosome","num.mark","left.index","right.index","seg.mean"),names(segments))
  if(!all(chk)) stop("Error in names of columns of data")
  if(!is.element(class(centromere),"data.frame")) stop("centromere info is not a data.frame")
  cchk<-is.element(c("chrom","left.base","right.base"),names(centromere))
  if(!all(cchk)) stop("Error in names of centromere data.frame")
  centromere$chrom[is.element(centromere$chrom, "X")] <- XchromCode(intenData)
  centromere$chrom[is.element(centromere$chrom, "Y")] <- YchromCode(intenData)
  centromere$chrom[is.element(centromere$chrom, "XY")] <- XYchromCode(intenData)
  centromere$chrom <- as.integer(centromere$chrom)

  # internal functions require these names, so convert from package standard
  names(segments)[names(segments) == "chromosome"] <- "chrom"
  names(segments)[names(segments) == "left.index"] <- "left"
  names(segments)[names(segments) == "right.index"] <- "right"

  ##delete any segments from male samples on chromosome X
  male<-sid[is.element(sex,"M")] 
  wdel<-is.element(segments$scanID,male)&segments$chrom==XchromCode(intenData)
  anoms<-segments[!wdel,]

  #find unsegmented chromosomes 
  smpchr<-paste(anoms$scanID,anoms$chrom)
  dup<-which(duplicated(smpchr))
  sng<-which(!is.element(smpchr,smpchr[dup]))
  anoms.sngl<-anoms[sng,]  #unsegmented chromosomes

  ######
  anoms2<-NULL  #raw with sd factor
  anoms.fil<-NULL #filtered
  normal.info<-NULL
  seg.info<-NULL
  samples<-unique(anoms$scanID)
  NS<-length(samples)
  for(i in 1:NS){
     snum<-samples[i]
     if(floor(i/10)*10-i==0 & verbose==TRUE){
       message(paste("processing ",i,"th scanID out of ",NS,sep=""))
     }

    sindex <- which(is.element(sid, snum))
    if(length(sindex)==0) stop(paste("Sample ",snum, " does not exist",sep=""))

    GENO <- getGenotype(genoData, snp=c(1,-1), scan=c(sindex,1))

    ## compute baf metric 
    sel<-is.element(intid,snp.ids) & (GENO == 1 | is.na(GENO))

    baf <- getBAlleleFreq(intenData, snp=c(1,-1), scan=c(sindex,1))
    ws<-!is.na(baf)

    INDEX<-which(sel&ws)
    CHR<-chrom[sel&ws]
    BAF<-baf[sel&ws]

    index<-NULL
    chr<-NULL
    baf.dat<-NULL
    uuch<-unique(anoms$chrom)
    uch<-uuch[uuch!=XYchromCode(intenData)]
    for(ch in uch){
      wc<-CHR==ch    #T/F for indices of CHR which match indices of INDEX,BAF
      bf<-BAF[wc]
      ind<-INDEX[wc]
      chrm<-CHR[wc]
      med<-median(bf,na.rm=TRUE)
      bf1<-1-bf
      bfm<-abs(bf-med)
      c<-cbind(bf,bf1,bfm)
      met<-apply(c,1,min)
      baf.metric<-sqrt(met)
      index<-c(index,ind)
      chr<-c(chr,chrm)
      baf.dat<-c(baf.dat,baf.metric)   
    } #end of chrom loop

    an<-anoms[is.element(anoms$scanID,snum),]
    an.sngl<-anoms.sngl[is.element(anoms.sngl$scanID,snum),]
  
    sel.chr.all<-an.sngl$chrom  
    if(sum(duplicated(an.sngl$chrom))!=0) stop(paste("Error in singletons for Sample ",snum,sep=" "))
    ## treat X chromosome and XY differently
    # not included in baseline but do need to be compared
    wX<-which(sel.chr.all==XchromCode(intenData))
    wps<-which(sel.chr.all==XYchromCode(intenData)) 
    if(length(wX)==0 & length(wps)==0){ sel.chr<-sel.chr.all} else {sel.chr<-sel.chr.all[-union(wX,wps)]}

    if(length(sel.chr)<2) {w.selec<-which(!is.element(chr,c(XchromCode(intenData),XYchromCode(intenData))))} else { 
  
      ##compare each autosome seg.mean with baseline based on other autosomes
      an.snglo<-an.sngl[order(an.sngl$seg.mean,decreasing=TRUE),]
      an.snglo<-an.snglo[is.element(an.snglo$chrom,sel.chr),]
      sel.chro<-an.snglo$chrom  #decreasing order of seg.mean so work with largest mean sequentially
      N<-length(sel.chro)  #sel.chro are autosome chroms unsegmented
      
      flag<-1
      j<-1
      while(flag==1){ 
        w.selec<-is.element(chr,sel.chro[(j+1):N])
        bbase<-baf.dat[w.selec]
        base.mean<-mean(bbase,na.rm=TRUE)
        base.sd<-sd(bbase,na.rm=TRUE)
        mean.chk<-an.snglo$seg.mean[is.element(an.snglo$chrom,sel.chro[j])]
        if(abs(mean.chk-base.mean)/base.sd >sd.long){
          if((N-j)==1) {flag<-0;keep<-N}
          j<-j+1
        }  else {keep<- j:N; flag<-0} 
      }#end while
      w.selec<-is.element(chr,sel.chro[keep])
    } #end of else - selection of unsegmented chroms to use
  
    ##baseline now based on autosome unsegmented not identified as whole chrom anoms

    bbase<-baf.dat[w.selec]
    base.mean<-mean(bbase,na.rm=TRUE)
    base.sd<-sd(bbase,na.rm=TRUE)
    sd.fac<-abs(an$seg.mean-base.mean)/base.sd
    
    if(length(sel.chr)<2){chr.ct<-0} else {chr.ct<-length(sel.chro[keep])}
    normi<-data.frame(snum,base.mean,base.sd,chr.ct)
    names(normi)<-c("scanID","base.mean","base.sd","chr.ct")
    normal.info<-rbind(normal.info,normi)

    an$sd.fac<-sd.fac
    an$sex<-sex[sindex]

    anoms2<-rbind(anoms2,an)
    ##an is segment/sd.fac info for the given sample
    ## base.mean and base.sd are for the given sample
    ##anoms2 now contains segment info along with sd.fac info - accumulating over samples


    ########## Merging and low percentage: sample/chrom ################
    an<-an[order(an$chrom,an$left),]
    an.seg.info<-NULL

    CHR<-unique(an$chrom)
    CHR<-CHR[CHR!=XYchromCode(intenData)]

    LRR <- getLogRRatio(intenData, snp=c(1,-1), scan=c(sindex,1))
    
    an3.fil<-NULL
    for(ch in CHR){
      anch<-an[an$chrom==ch,]
      if(dim(anch)[1]==0) next
      nsegs<-dim(anch)[1]
      tp<-data.frame(snum,ch,nsegs)
      names(tp)<-c("scanID","chrom","num.segs")
      an.seg.info<-rbind(an.seg.info,tp)  

      tmp2<-GWASTools:::mergeSeg(anch,snum,ch,centromere,Pos,base.mean,base.sd,sd.reg,sd.long,num.mark.thresh,long.num.mark.thresh,low.frac.used) ## merged for given samp/chrom

      ### modifying breakpoints for potential homo deletions #####       
        ### delHomoRuns returns original breakpoints for any anom not needing adjustment

      tst.rev<-GWASTools:::delHomoRuns(tmp2,sid,snp.ids,intid,LRR,
        run.size,inter.size,low.frac.used,lrr.cut,ct.thresh,frac.thresh)

      diff.right<-tst.rev$old.right-tst.rev$right
      diff.left<-tst.rev$old.left-tst.rev$left
      d.keep<-NULL
      tst.rev$homodel.adjust<-FALSE
      for(k in 1:dim(tst.rev)[1]) {
        if(diff.right[k]==0 & diff.left[k]==0) next
        d.keep<-c(d.keep,k)
        int2<-intid>=intid[tst.rev$left[k]] & intid<=intid[tst.rev$right[k]]
        selgood<-is.element(intid,snp.ids)
        whm<- GENO == 1 | is.na(GENO)
        index2<-which(int2&selgood&whm)  
        mt<-is.element(index,index2)
        seg.mn<-mean(baf.dat[mt],na.rm=TRUE)
        sdf<-abs(seg.mn-base.mean)/base.sd
        tst.rev$num.mark[k]<-length(index2)
        tst.rev$sd.fac[k]<-sdf
        tst.rev$seg.mean[k]<-seg.mn  
      }
      if(length(d.keep)!=0) tst.rev$homodel.adjust[d.keep]<-TRUE
       # want to include small homo dels originally found exactly (so not adjusted)
      ct<-tst.rev$low.ct
      nlrr<-tst.rev$nmark.lrr
      ws<-which(ct>=ct.thresh & ct/nlrr>=.9) 
      #keeps homo dels found exactly (didn't need adjustment) that might be too small to pass further filters
      d.keep<-union(d.keep,ws)
      colm<-c("scanID","chrom","left","right","num.mark","seg.mean","sd.fac","sex","merge","homodel.adjust")
      tst.rev<-tst.rev[,colm]
     
 
      ######## FILTER ######################
      s1<-d.keep #will keep any homo dels found
      s2<-tst.rev$num.mark>=num.mark.thresh& tst.rev$sd.fac>=sd.reg
      s3<-tst.rev$num.mark>long.num.mark.thresh& tst.rev$sd.fac>=sd.long
      s<-which(s2|s3)
      filin<-union(s,d.keep)
      fil<-tst.rev[filin,]
      an3.fil<-rbind(an3.fil,fil) #becomes filtered anoms for current sample
    }#end ch loop

    anoms.fil<-rbind(anoms.fil,an3.fil)
    seg.info<-rbind(seg.info,an.seg.info)
  } #end of sample loop

  anoms2<-anoms2[order(anoms2$scanID,anoms2$chrom,anoms2$left),] #raw annotated
  anoms2$left.base <- getPosition(intenData, index=anoms2$left)
    anoms2$right.base <- getPosition(intenData, index=anoms2$right)

  if(!is.null(anoms.fil)){
    anoms.fil<-anoms.fil[order(anoms.fil$scanID,anoms.fil$chrom,anoms.fil$left),]

    # convert index to base position
    anoms.fil$left.base <- getPosition(intenData, index=anoms.fil$left)
    anoms.fil$right.base <- getPosition(intenData, index=anoms.fil$right)
  }

  ## XY filter
  tmp<-anoms2[anoms2$chrom==XYchromCode(intenData),]
  if(dim(tmp)[1]!=0){
    s1<-tmp$num.mark>=num.mark.thresh& tmp$sd.fac>=sd.reg
    s2<-tmp$num.mark>long.num.mark.thresh& tmp$sd.fac>=sd.long
    if(sum(s1|s2)!=0){
      an.XY<-tmp[s1|s2,]
      an.XY$merge<-NA
      an.XY$homodel.adjust<-NA
      an.XY$left.base<-getPosition(intenData, index=an.XY$left)
      an.XY$right.base<-getPosition(intenData, index=an.XY$right)} else an.XY<-NULL
    anoms.fil<-rbind(anoms.fil,an.XY)
    an.seg.info<-NULL
    s<-unique(tmp$scanID)
    for(snum in s){ 
      an<-tmp[tmp$scanID==snum,]
      ns<-dim(an)[1]
      tt<-data.frame(snum,XYchromCode(intenData),ns)
      names(tt)<-c("scanID","chrom","num.segs")
      an.seg.info<-rbind(an.seg.info,tt)
    }
    seg.info<-rbind(seg.info,an.seg.info)
  }
 
  ## further filtering based on low.frac.used and/or messiness
  if(!is.null(anoms.fil)& dim(anoms.fil)[1]!=0){
    nf<-dim(anoms.fil)[1]  
    frac.used<-NULL
    for(kk in 1:nf){
      int<-intid>=intid[anoms.fil$left[kk]] & intid<=intid[anoms.fil$right[kk]]
      selgood<-is.element(intid,snp.ids)
      int.ok<-int&selgood
      frac.used<-c(frac.used,anoms.fil$num.mark[kk]/sum(int.ok)) 
    }
    anoms.fil$frac.used<-frac.used 
    wd1<-anoms.fil$frac.used<=very.low.frac.used
    wd2<-anoms.fil$frac.used<low.frac.used & anoms.fil$num.mark<low.frac.used.num.mark
    if(!is.null(low.qual.ids)){
      wd3<-anoms.fil$frac.used<low.frac.used & anoms.fil$num.mark<low.qual.frac.num.mark & is.element(anoms.fil$scanID,low.qual.ids)
    }
    wdel<-wd1 | wd2 
    if(!is.null(low.qual.ids)) wdel<-wdel|wd3
    anoms.fil<-anoms.fil[!wdel,]
  }
  seg.info<-seg.info[order(seg.info$scanID,seg.info$chrom),]
  out<-list(anoms2,anoms.fil,normal.info,seg.info)
  names(out)<-c("raw","filtered","base.info","seg.info")
  
  # convert back to package standard names
  for (i in 1:length(out)) {
    if ("chrom" %in% names(out[[i]]))
      names(out[[i]])[names(out[[i]]) == "chrom"] <- "chromosome"
    if ("left" %in% names(out[[i]]))
      names(out[[i]])[names(out[[i]]) == "left"] <- "left.index"
    if ("right" %in% names(out[[i]]))
      names(out[[i]])[names(out[[i]]) == "right"] <- "right.index"
  }

  return(out)
} #end of function


