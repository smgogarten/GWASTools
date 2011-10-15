#####
# Start from stored baf/lrr sd data
# Exclude low quality and no.post samples from baf/lrr sd 
# Return the medians of autosomal baf/lrr sd, one per sample
#####

## Revised by Cecelia Laurie to include information about sample number

medianSdOverAutosomes <- function(sd.by.scan.chrom.window)
{
  isd <- sd.by.scan.chrom.window

  snums<-as.integer(rownames(isd[[1]]))
  # include autosomes only
  isd.df.sub <- isd[names(isd) %in% 1:22] #
  # get median sd for autosomes. one value per sample
  isd.df.sub <- data.frame(isd.df.sub, stringsAsFactors=FALSE)

  med.sd.per.samp <- apply(isd.df.sub, 1, median)
  med.sd.samp.info<-data.frame(snums,med.sd.per.samp)
  names(med.sd.samp.info)<-c("scanID","med.sd")
  return(med.sd.samp.info)
  
}
