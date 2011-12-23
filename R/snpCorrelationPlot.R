
snpCorrelationPlot<-function(correlations,
                               chromosome,
                               chrom.labels=c(1:22,"X","XY","Y","M"),
                               ylim=c(0,1),
                               ylab="abs(correlation)",
                               ...){

	# correlations is a vector of correlations between SNPs and a given PC
	# chromosome is a parallel vector of chromosome IDs (chrom.int)
	# v2 makes pch="." and removes cex=0.3

	nc <- length(unique(chromosome))
	# in case the chromosomes are not all  consecutive integers
	chromosome <- as.numeric(as.factor(chromosome))
		
	if(nc!=length(chrom.labels)) stop("chrom.labels length not equal to number of chromosomes")
	N <- length(correlations)

        chromstart <- which(c(1,diff(chromosome)) == 1) # element numbers for the first SNP on each chromosome
        chromend <- c(chromstart[-1],N)  # element numbers for the last SNP on each chromosome

        x <- (1:N) + chromosome*(chromend[1]/6) # uniform spacing for SNP positions with offset
	y <- correlations
        
        plot(x, y, col=chromosome, pch=".",
             bg=chromosome, xlab="Chromosome", ylab=ylab,
             xaxt="n", ..., ylim=ylim)

        centers <- (x[chromstart]+x[chromend]-(chromend[1]/6))/2
	centers[length(chrom.labels)] <- x[N]
        axis(1, at=centers, label=chrom.labels, las=2)
    }



