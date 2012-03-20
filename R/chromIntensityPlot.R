chromIntensityPlot <- function(
		 intenData, # object of type IntensityData
                 scan.ids, 
		 chrom.ids, 
		 type = c("BAF/LRR", "BAF", "LRR", "R", "Theta", "R/Theta"),
		 code = NULL, 
		 main.txt = NULL, 
		 abln = NULL, 
		 horizln = c(1/2, 1/3, 2/3), 
		 colorGenotypes = FALSE, 
		 genoData = NULL,
                 colorBatch = FALSE,
                 batch.column = NULL,
		 snp.exclude = NULL,
                 cex=0.5,
		 ...) 
{
    # check arguments
    if (length(scan.ids) != length(chrom.ids)) {
        stop("scan.ids and chrom.ids must be parallel vectors of the same length")
      }
    type <- match.arg(type)
    if (type == "BAF" | type == "BAF/LRR") {
        if (!hasBAlleleFreq(intenData)) stop("BAlleleFreq not found")
      }
    if (type == "LRR" | type == "BAF/LRR") {
        if (!hasLogRRatio(intenData)) stop("LogRRatio not found")
      }
    if (type == "R" | type == "Theta" | type == "R/Theta") {
        if (!hasX(intenData) | !hasY(intenData)) stop("X and Y not found")
      }
    if(!is.null(code)) { stopifnot(length(code)==length(scan.ids)) }
    if (colorGenotypes & colorBatch) {
        stop("cannot color by Genotype and Batch simultaneously")
      }
    if (colorGenotypes & is.null(genoData)) {
        stop("genoData must be specified if colorGenotypes is TRUE")
    }
    if (colorBatch) {
      if (is.null(batch.column)) {
        stop("batch.column must be specified if colorBatch is TRUE")
      } else if (!hasSnpVariable(intenData, batch.column)) {
        stop(paste("SNP variable", batch.column, "was not found in intenData."))
      }
    }
    
    chr <- getChromosome(intenData)
    chr.char <- getChromosome(intenData, char=TRUE)
    pos <- getPosition(intenData)
    scanID <- getScanID(intenData)
    
    for (i in 1:length(scan.ids)) {
	# index for sample to be plotted
        stopifnot(is.element(scan.ids[i], scanID))
        sid <- which(scanID == scan.ids[i])
        
        chri <- which(is.element(chr, chrom.ids[i]))
        chr.start <- chri[1]
        chr.count <- length(chri)
        chr.label <- unique(chr.char[chri])
        
	# logical for snps to be plotted
        if (!is.null(snp.exclude)) {
            snpID <- getSnpID(intenData, index=chri)
            toPlot <- !(snpID %in% snp.exclude)
        } else {
            toPlot <- TRUE
        }
        
        # get the data to plot
        if (type == "BAF" | type == "LRR" | type == "BAF/LRR") {
            logrratio <- getLogRRatio(intenData, snp=c(chr.start,chr.count), scan=c(sid,1))
            mninten <- mean(logrratio, na.rm = T)
            bafreq <- getBAlleleFreq(intenData, snp=c(chr.start,chr.count), scan=c(sid,1))
        } else { # type is R, Theta or R/Theta
            x <- getX(intenData, snp=c(chr.start,chr.count), scan=c(sid,1))
            y <- getY(intenData, snp=c(chr.start,chr.count), scan=c(sid,1))
            theta <- atan(y/x) * (2/pi)
            r <- x + y
        }
        posi <- pos[chri]
    
	# calculate values at which the 1/8 lines should be plotted
        top <- max(posi, na.rm = TRUE)
        bot <- min(posi, na.rm = TRUE)
        len <- top - bot
        eighth <- len/8
        c <- 1:8
        vals <- rep(eighth * c) + bot

        # plot title
        if (is.null(main.txt)) {
            txt <- paste("Scan", scan.ids[i], "- Chromosome", chr.label)
            if (!is.null(code)) {
              txt <- paste(txt, "-", code[i])
            }
        } else {
            txt <- main.txt
        }

        # create genotype color vector if colorGenotypes==TRUE
        gcol <- rep("black", length(chri))
        if (colorGenotypes) {
            genos <- getGenotype(genoData, snp=c(chr.start,chr.count), scan=c(sid,1))
            gcol[!is.na(genos) & genos == 0] <- "blue"
            gcol[!is.na(genos) & genos == 1] <- "green"
            gcol[!is.na(genos) & genos == 2] <- "red"
        }

        # create batch color vector if colorBatch==TRUE
        if (colorBatch) {
          batch <- getSnpVariable(intenData, batch.column, index=chri)
          uniqBatch <- unique(batch)
          numpool <- length(uniqBatch)
          col <- colors()[seq(from=450,to=637,length.out=numpool)]
          for(i in 1:length(col)) gcol[is.element(batch,uniqBatch[i])] <- col[i]
        }
    
	# make the plots
        if (!is.null(abln)) {
            abst <- abln[(i * 2) - 1]
            aben <- abln[i * 2]
        }
        else {
            abst <- -1
            aben <- -1
        }
        if (!is.null(horizln)) {
            hv <- vector()
            horizvals <- format(horizln, digits = 4)
            for (g in 1:length(horizvals)) {
                hv <- paste(hv, horizvals[g])
            }
            subnm <- paste("horizontal line =", hv)
        }
        if (type == "BAF/LRR" | type=="R/Theta") {
            par(mfrow=c(2,1))
        }
        if (type == "LRR" | type == "BAF/LRR") {
            plot((posi/1e+06)[toPlot], logrratio[toPlot], xlab = "position (Mb)", 
                ylab = "LRR", sub = "horizontal line = mean LRR", 
                main = txt, ylim = c(-2, 2),
                type = "n", ...)
            for (d in 1:length(vals)) {
                abline(v = vals[d]/1e+06, col = "royalblue", 
                  lty = 3, lwd = 2)
            }
            if (abst != -1) 
                abline(v = abst, col = "red", lty = 2, lwd = 1.2)
            if (aben != -1) 
                abline(v = aben, col = "red", lty = 2, lwd = 1.2)
            points((posi/1e+06)[toPlot], logrratio[toPlot], cex=cex, ...)
            abline(h = mninten, col = "gray")
        }
        if (type == "BAF" | type == "BAF/LRR") {
            plot((posi/1e+06)[toPlot], bafreq[toPlot], type = "n", 
                xlab = "position (Mb)", ylab = "BAF", 
                sub = subnm, main = txt, col = gcol[toPlot], ...)
            if (abst != -1) 
                abline(v = abst, col = "red", lty = 2, lwd = 1.3)
            if (aben != -1) 
                abline(v = aben, col = "red", lty = 2, lwd = 1.3)
            for (d in 1:length(vals)) {
                abline(v = vals[d]/1e+06, col = "royalblue", 
                  lty = 3, lwd = 2)
            }
            points((posi/1e+06)[toPlot], bafreq[toPlot], col = gcol[toPlot], cex=cex, ...)
            abline(h = horizln, col = "gray")
        }
        if (type == "R" | type == "R/Theta") {
            plot((posi/1e+06)[toPlot], r[toPlot], xlab = "position (Mb)", 
                ylab = "R", main = txt, type = "n", ...)
            for (d in 1:length(vals)) {
                abline(v = vals[d]/1e+06, col = "royalblue", 
                  lty = 3, lwd = 2)
            }
            if (abst != -1) 
                abline(v = abst, col = "red", lty = 2, lwd = 1.3)
            if (aben != -1) 
                abline(v = aben, col = "red", lty = 2, lwd = 1.3)
            points((posi/1e+06)[toPlot], r[toPlot], col = gcol[toPlot], cex=cex, ...)
        }
        if (type == "Theta" | type == "R/Theta") {
            plot((posi/1e+06)[toPlot], theta[toPlot], xlab = "position (Mb)", 
                ylab = "Theta", main = txt, ...)
            for (d in 1:length(vals)) {
                abline(v = vals[d]/1e+06, col = "royalblue", 
                  lty = 3, lwd = 2)
            }
            if (abst != -1) 
                abline(v = abst, col = "red", lty = 2, lwd = 1.3)
            if (aben != -1) 
                abline(v = aben, col = "red", lty = 2, lwd = 1.3)
            points((posi/1e+06)[toPlot], r[toPlot], col = gcol[toPlot], cex=cex, ...)
        }
    }
}
