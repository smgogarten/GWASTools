manhattanPlot <- function(p,
                          chromosome,
                          ylim = NULL,
                          trunc.lines = TRUE,
                          signif = 5e-8,
                          thinThreshold=NULL,
                          pointsPerBin=10000,
                          col=NULL,
                          ...)
{
  stopifnot(length(p) == length(chromosome))
  if (!is.null(col)) stopifnot(length(p) == length(col))
  
  logp <- -log(p,10)
  
  # remove NAs
  sel <- !is.na(logp)
  logp <- logp[sel]
  chromosome <- chromosome[sel]
  if (!is.null(col)) col <- col[sel]
  
  # select SNPs if thinning is requested
  if (!is.null(thinThreshold)){
    #quant <- quantile(logp[logp < thinThreshold], probs=1:10/10)
    breaks <- seq(from=0, to=thinThreshold, length.out=11)
    logp.cut <- cut(logp, breaks=breaks, right=FALSE)
    ind.list <- list()
    for (level in levels(logp.cut)){
      sel <- which(logp.cut == level)
      ind.list[[level]] <- sel[sample.int(length(sel), min(1000, length(sel)))]
    }
    ind.list[["max"]] <- which(logp >= thinThreshold)
    ind <- unlist(ind.list, use.names=FALSE)
    ind <- sort(ind) # sorting necessary for chromosome ordering
    
    p <- p[ind]
    logp <- logp[ind]
    chromosome <- chromosome[ind]
    if (!is.null(col)) col <- col[ind]
  } 
  
  N <- length(logp)
  ymax <- log(N,10) + 4
  if (is.null(ylim)) {
    ylim <- c(0,ymax)
  } else { 
    ymax <- ylim[2]
  }
  
  chrom.labels <- unique(chromosome)
  chrom.int <- as.numeric(factor(chromosome, levels=chrom.labels))
  chromstart <- which(c(1,diff(chrom.int)) == 1) # element numbers for the first SNP on each chromosome
  chromend <- c(chromstart[-1],N)  # element numbers for the last SNP on each chromosome
  x <- (1:N) + chrom.int*(chromend[1]/6) # uniform spacing for SNP positions with offset
  
  if (is.null(col)){
  # colors from colorbrewer Dark2 map, 8 levels
    chrom.col <- factor(chrom.int)
    levels(chrom.col) <- rep_len(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), length(levels(chrom.col)))
    col <- as.character(chrom.col)
  }
  
  trunc <- FALSE
  if (trunc.lines & any(logp > ymax, na.rm=TRUE)) trunc <- TRUE
  y <- pmin(ymax,logp)
  
  plot(x, y, cex=0.3+y/ymax, col=col, pch=ifelse(y==ymax,24,19),
       bg=chrom.int, xlab="Chromosome", ylab=quote(-log[10](p)),
       xaxt="n", ..., ylim=ylim)
  if (trunc) lines(c(min(x), max(x)), c(ymax,ymax), col = 'grey')
  
  centers <- (x[chromstart]+x[chromend]-(chromend[1]/6))/2
  centers[length(chrom.labels)] <- x[N]  # puts the last tick at the position of last snp
  axis(1, at=centers, labels=chrom.labels, las=2)
  
  # add a line for genome-wide significance
  if (!is.null(signif)) {
    abline(h=-log10(signif), lty=2, col="gray")
  }
}

