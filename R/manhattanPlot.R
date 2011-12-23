
manhattanPlot <- function(p,
		 chromosome,
		 chrom.labels=c(1:22,"X","XY","Y","M"),
		 ylim = NULL,
		 trunc.lines = TRUE,
		 ...)
{

	nc <- length(unique(chromosome))
	# in case the chromosomes are not all  consecutive integers
	chromosome <- as.numeric(as.factor(chromosome))
		
	if(nc!=length(chrom.labels)) stop("chrom.labels length not equal to number of chromosomes")
   	logp <- -log(p,10)

    	N <- length(logp)
   	ymax <- log(N,10) + 4
   	if (is.null(ylim)) {
          ylim <- c(0,ymax)
        } else { 
          ymax <- ylim[2]
        }
   	
   	chromstart <- which(c(1,diff(chromosome))==1) # element numbers for the first SNP on each chromosome
        chromend <- c(chromstart[-1],N)  # element numbers for the last SNP on each chromosome

   	x <- (1:N) + chromosome*(chromend[1]/6) # uniform spacing for SNP positions with offset
   	trunc <- FALSE
   	if (trunc.lines & any(logp > ymax, na.rm=TRUE)) trunc <- TRUE
    	y <- pmin(ymax,logp)

        plot(x, y, cex=0.3+y/ymax, col=chromosome, pch=ifelse(y==ymax,24,19),
             bg=chromosome, xlab="Chromosome", ylab=quote(-log[10](p)),
             xaxt="n", ..., ylim=ylim)
        if (trunc) lines(c(min(x), max(x)), c(ymax,ymax), col = 'grey')

        centers <- (x[chromstart]+x[chromend]-(chromend[1]/6))/2
        centers[nc] <- x[N]  # puts the last tick at the position of last snp
        axis(1, at=centers, label=chrom.labels, las=2)
	
 }


