# Description: This function determines which scans/chromosomes
# are a given number of SDs from the mean, using new BAF values
# produced by meanBAFSDbyChromWindow. The function returns a matrix
# (with columns "scanID", "chromosome", "bin", and "sex")
# of scan/chromosome combinations more than X SDs from the mean.

findBAFvariance <- function(sd.by.chrom.window, sd.by.scan.chrom.window,
                            sex, sd.threshold)
{
  chromosomes <- names(sd.by.scan.chrom.window)
  n.chromosomes <- length(chromosomes)  
  
  # create a matrix to hold scan, chr combination that satisfies our SDs from mean comparison
  res <- vector()

  r <- 0

  # loop through chromosomes
  for(s in 1:n.chromosomes) {
    # autosomes
    if (chromosomes[s] != "X") {
      # loop through scans
      for(n in 1:nrow(sd.by.scan.chrom.window[[s]])) { 
        bc <- 0
        # loop through bins
        for(b in 1:ncol(sd.by.scan.chrom.window[[s]])) { 
          sdev <- sd.by.chrom.window[[s]]["SD",b]
          m <- sd.by.chrom.window[[s]]["Mean",b]
          val <- sd.by.scan.chrom.window[[s]][n,b]
          samp <- rownames(sd.by.scan.chrom.window[[s]])[n]
          if(!is.na(val) & !is.na(sdev) & !is.na(m) &
             val > sd.threshold*sdev+m) {
            bc <- bc+1
            if(bc==1) { 
              res <- rbind(res, c(samp, chromosomes[s], bc, sex[n]))
              r <- r+1
            } else {
              res[r,3] <- bc
            } 
          }
        } # end loop through bins
      } # end loop through scans
    }  # end if autosome

    # X chrom
    if (chromosomes[s] == "X") {
      # loop through scans
      for(n in 1:nrow(sd.by.scan.chrom.window[[s]])) { 
        bc <- 0
        for(b in 1:ncol(sd.by.scan.chrom.window[[s]])) {
          # loop through bins
          sdf <- sd.by.chrom.window[[s]]["Female SD",b]
          mf <- sd.by.chrom.window[[s]]["Female Mean",b]
          sdm <- sd.by.chrom.window[[s]]["Male SD",b]
          mm <- sd.by.chrom.window[[s]]["Male Mean",b]
          val <- sd.by.scan.chrom.window[[s]][n,b]
          samp <- rownames(sd.by.scan.chrom.window[[s]])[n]
          if(sex[n] == "F") { # it's a female
            if(!is.na(val) & !is.na(sdf) & !is.na(mf) &
               val > sd.threshold*sdf+mf) {
              bc <- bc+1
              if(bc==1) { 
                res <- rbind(res, c(samp, "X", bc, sex[n]))
                r <- r+1
              } else {
                res[r,3] <- bc
              } 
            }
          } else { # it's a male
            if(!is.na(val) & !is.na(sdm) & !is.na(mm) &
               val > sd.threshold*sdm+mm) { 
              bc <- bc+1
              if(bc==1) { 
                res <- rbind(res, c(samp, "X", bc, sex[n]))
                r <- r+1
              } else {
                res[r,3] <- bc
              } 
            } 
          }
	} # end loop through bins
      } # end loop through scans
    } # end if X
  } # end loop through chromosomes
  colnames(res) <- c("scanID", "chromosome", "bin", "sex")
  return(res)
}
