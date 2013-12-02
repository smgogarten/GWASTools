
snpCorrelationPlot<-function(correlations,
                             chromosome,
                             ylim=c(0,1),
                             ylab="abs(correlation)",
                             ...){
  
  # correlations is a vector of correlations between SNPs and a given PC
  # chromosome is a parallel vector of chromosome IDs (chrom.int)
  # v2 makes pch="." and removes cex=0.3
  
  stopifnot(length(correlations) == length(chromosome))
  
  N <- length(correlations)
  y <- correlations
  
  chrom.labels <- unique(chromosome)
  chrom.int <- as.numeric(factor(chromosome, levels=chrom.labels))
  chromstart <- which(c(1,diff(chrom.int)) == 1) # element numbers for the first SNP on each chromosome
  chromend <- c(chromstart[-1],N)  # element numbers for the last SNP on each chromosome
  x <- (1:N) + chrom.int*(chromend[1]/6) # uniform spacing for SNP positions with offset
  
  # colors from colorbrewer Dark2 map, 8 levels
  chrom.col <- factor(chrom.int)
  levels(chrom.col) <- rep_len(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), length(levels(chrom.col)))
  chrom.col <- as.character(chrom.col)
  
  
  plot(x, y, col=chrom.col, pch=".",
       bg=chrom.int, xlab="Chromosome", ylab=ylab,
       xaxt="n", ..., ylim=ylim)
  
  centers <- (x[chromstart]+x[chromend]-(chromend[1]/6))/2
  centers[length(chrom.labels)] <- x[N] # puts the last tick at the position of last snp
  axis(1, at=centers, labels=chrom.labels, las=2)
}



