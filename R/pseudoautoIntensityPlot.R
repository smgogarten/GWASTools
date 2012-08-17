#####


pseudoautoIntensityPlot <- function(intenData, # object of type IntensityData
                                    scan.ids,
                                    main = NULL,
                                    plotY = FALSE,
                                    hg.build = c("hg18", "hg19"),
                                    snp.exclude = NULL,
                                    cex=0.5,
                                    ...)
{
# plot BAF/LRR plots for sample nums listed in scan.ids, overlaying the X and XY snps
  
  if (!hasLogRRatio(intenData) | !hasBAlleleFreq(intenData)) stop("required variables not found")
  
  if (!is.null(main)) {    
    if (length(main) == 1 & length(scan.ids) > 1) {
      main <- rep(main, length(scan.ids))
    } else {
      stopifnot(length(main) == length(scan.ids))
    }
  }

  scanID <- getScanID(intenData)
  hasSexLabel <- hasSex(intenData)
  if (hasSexLabel) sex <- getSex(intenData)

  chr <- getChromosome(intenData, char=TRUE)
  pos <- getPosition(intenData)
  sexChrom <- which(chr %in% c("X", "Y", "XY"))
  chr.start <- sexChrom[1]
  chr.count <- length(sexChrom)
  chr <- chr[sexChrom]
  pos <- pos[sexChrom]

  logrratio <- matrix(nrow=chr.count, ncol=length(scan.ids))
  bafreq <- matrix(nrow=chr.count, ncol=length(scan.ids))
  sexLabel <- rep("", length(scan.ids))
  for(i in 1:length(scan.ids)) {
    thisi <- which(scanID == scan.ids[i])
    if (hasSexLabel) sexLabel[i] <- paste("Sex", sex[thisi])
    logrratio[,i] <- getLogRRatio(intenData, snp=c(chr.start, chr.count), scan=c(thisi,1))
    bafreq[,i] <- getBAlleleFreq(intenData, snp=c(chr.start, chr.count), scan=c(thisi,1))
  }

  if (!is.null(snp.exclude)) {
    snpID <- getSnpID(intenData, index=sexChrom)
    keep <- !(snpID %in% snp.exclude)
    logrratio <- logrratio[keep,,drop=FALSE]
    bafreq <- bafreq[keep,,drop=FALSE]
    chr <- chr[keep]
    pos <- pos[keep]
  }
  
  # get the X, Y and XY inds
  xinds <- !is.na(chr) & chr=="X"
  xyinds <- !is.na(chr) & chr=="XY"
  yinds <- !is.na(chr) & chr=="Y"
  xcol <- rep(NA,length(chr))
  xcol[xinds] <- "magenta"
  xcol[yinds] <- "skyblue"
  xcol[xyinds] <- "darkgreen"

  allinds <- xinds+xyinds+yinds
  allinds <- allinds==1

  xxyinds <- xinds+xyinds
  xxyinds <- xxyinds==1 # so it's a logical vector
  yxyinds <- yinds+xyinds
  yxyinds <- yxyinds==1 # so it's a logical vector
  xcolxxy <- xcol[xxyinds]
  xcolyxy <- xcol[yxyinds]

  # rectangles with XY regions (in Mb)
  hg.build <- match.arg(hg.build)
  pa <- get(data(list=paste("pseudoautosomal", hg.build, sep="."),
                 package="GWASTools", envir=environment()))
  PAR1start <- pa["X.PAR1", "start.base"] / 1e6
  PAR1end <- pa["X.PAR1", "end.base"] / 1e6
  xXTRstart <-  pa["X.XTR", "start.base"] / 1e6
  xXTRend <- pa["X.XTR", "end.base"] / 1e6
  yXTRstart <- pa["Y.XTR", "start.base"] / 1e6
  yXTRend <- pa["Y.XTR", "end.base"] / 1e6
  xPAR2start <-  pa["X.PAR2", "start.base"] / 1e6
  xPAR2end <- pa["X.PAR2", "end.base"] / 1e6
  yPAR2start <- pa["Y.PAR2", "start.base"] / 1e6
  yPAR2end <- pa["Y.PAR2", "end.base"] / 1e6

  # set some plotting parameters
  if (plotY) {
    par(mfcol=c(2,2))
  } else {
    par(mfrow=c(2,1))
  }

  for(i in 1:length(scan.ids)) {

        # X chromosome
        if(!is.null(main)) {
          txt.main <- main[i]
        } else {
          txt.main <- paste("Scan",scan.ids[i],"-",sexLabel[i])
        }
        txt.leg <- paste("magenta = X SNPs, green = XY SNPs",
                         "gray = PAR1/PAR2, yellow=XTR", sep="\n")

        # LRR
        plot(pos[xxyinds]/1e6,logrratio[xxyinds,i],xlab="position (Mb)",ylab="LRR",sub="horizontal line = mean LRR",main=txt.main,col=xcolxxy,type="n", ylim=c(-2,2), ...)
        mninten <- mean(logrratio[xxyinds,i],na.rm=TRUE)
        abline(h=mninten,col="gray")

        # expected XY rectangles
        rect(xleft=c(PAR1start, xPAR2start), xright=c(PAR1end, xPAR2end),
             ybottom=c(-2,-2), ytop=c(2,2), col="gray", border=NA)
        rect(xleft=xXTRstart, xright=xXTRend,
             ybottom=c(-2,-2), ytop=c(2,2), col="yellow", border=NA)
        
        # overlay the xy points on the plot
        points(pos[xinds]/1e6,logrratio[xinds,i],col="magenta", cex=cex)
        points(pos[xyinds]/1e6,logrratio[xyinds,i],col="darkgreen", cex=cex)

        # BAF
        plot(pos[xxyinds]/1e6,bafreq[xxyinds,i],xlab="position (Mb)",ylab="BAF",main=txt.leg,col=xcolxxy,type="n", ...)

        # expected XY rectangles
        rect(xleft=c(PAR1start, xPAR2start), xright=c(PAR1end, xPAR2end),
             ybottom=c(0,0), ytop=c(1,1), col="gray", border=NA)
        rect(xleft=xXTRstart, xright=xXTRend,
             ybottom=0, ytop=1, col="yellow", border=NA)
        
        # overlay the xy points on the plot
        points(pos[xinds]/1e6,bafreq[xinds,i],col="magenta", cex=cex)
        points(pos[xyinds]/1e6,bafreq[xyinds,i],col="darkgreen", cex=cex)

        # Y chromosome
        if (plotY) {
          txt.leg <- paste("blue=Y SNPs, green=XY SNPs",
                           "gray = PAR1/PAR2, yellow=XTR", sep="\n")

          # scale base positions for XY regions to Y chrom
          ypts <- pos[yinds]/1e6
          ylrr <- logrratio[yinds,i]
          ybaf <- bafreq[yinds,i]
          
          xypts <- pos[xyinds]/1e6
          xylrr <- logrratio[xyinds,i]
          xybaf <- bafreq[xyinds,i]

          xtrsel <- xypts >= (xXTRstart - 1) & xypts <= (xXTRend + 1)
          xypts[xtrsel] <- xypts[xtrsel] - (xXTRstart - yXTRstart)

          par2sel <- xypts >= xPAR2start
          xypts[par2sel] <- xypts[par2sel] - (xPAR2end - yPAR2end)

  
          # LRR
          plot(c(ypts, xypts), c(ylrr, xylrr),xlab="position (Mb)",ylab="LRR",sub="horizontal line = mean LRR",main=txt.main,col=xcolyxy,type="n", ylim=c(-2,2), ...)
          mninten <- mean(logrratio[yxyinds,i],na.rm=TRUE)
          abline(h=mninten,col="gray")

          # expected XY rectangles
          rect(xleft=c(PAR1start, yPAR2start), xright=c(PAR1end, yPAR2end),
               ybottom=c(-2,-2), ytop=c(2,2), col="gray", border=NA)
          rect(xleft=yXTRstart, xright=yXTRend,
               ybottom=c(-2,-2), ytop=c(2,2), col="yellow", border=NA)
        
          # overlay the xy points on the plot
          points(ypts, ylrr, col="skyblue", cex=cex)
          points(xypts, xylrr,col="darkgreen", cex=cex)

          # BAF
          plot(c(ypts, xypts), c(ybaf, xybaf),xlab="position (Mb)",ylab="BAF",main=txt.leg,col=xcolyxy,type="n", ...)

          # expected XY rectangles
          rect(xleft=c(PAR1start, yPAR2start), xright=c(PAR1end, yPAR2end),
               ybottom=c(0,0), ytop=c(1,1), col="gray", border=NA)
          rect(xleft=yXTRstart, xright=yXTRend,
               ybottom=0, ytop=1, col="yellow", border=NA)
        
          # overlay the xy points on the plot
          points(ypts, ybaf, col="skyblue", cex=cex)
          points(xypts, xybaf, col="darkgreen", cex=cex)
        }
  }

}
