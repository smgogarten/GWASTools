chromIntensityPlot <- function(
		 intenData, # object of type IntensityData
                 scan.ids,
		 chrom.ids,
		 type = c("BAF/LRR", "BAF", "LRR", "R", "Theta", "R/Theta"),
		 main = NULL,
		 info = NULL,
		 abln = NULL,
		 horizln = c(1/2, 1/3, 2/3),
		 colorGenotypes = FALSE,
		 genoData = NULL,
                 colorBatch = FALSE,
                 batch.column = NULL,
		 snp.exclude = NULL,
                 ideogram=TRUE, ideo.zoom=TRUE, ideo.rect=FALSE,
                 cex=0.5,
                 cex.leg=1.5,
                 colors = c("default", "neon", "primary"),
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
    if (!is.null(main)) {
      if (length(main) == 1 & length(scan.ids) > 1) {
        main <- rep(main, length(scan.ids))
      } else {
        stopifnot(length(main) == length(scan.ids))
      }
    }
    if (!is.null(info)) { stopifnot(length(info)==length(scan.ids)) }
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
    colors <- match.arg(colors)

    chr <- getChromosome(intenData)
    chr.char <- getChromosome(intenData, char=TRUE)
    pos <- getPosition(intenData)
    scanID <- getScanID(intenData)
    if (hasSex(intenData)) {
      sex <- getSex(intenData)
    } else {
      sex <- NULL
    }

    chrom.char <- chrom.ids
    chrom.char[chrom.ids == XchromCode(intenData)] <- "X"
    chrom.char[chrom.ids == YchromCode(intenData)] <- "Y"

    # layout
    if (ideogram & type %in% c("BAF/LRR", "R/Theta")) {
      layout(matrix(c(1,2,3), nrow=3, ncol=1), heights=c(0.4, 0.4, 0.2))
    } else if (ideogram | type %in% c("BAF/LRR", "R/Theta")) {
      layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(0.5, 0.5))
    }

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
            toPlot <- rep(TRUE, length(chri))
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
        if (is.null(main)) {
            if (!is.null(sex)) {
              txt.main <- paste("Scan", scan.ids[i], "-", sex[sid], "- Chromosome", chr.label)
            } else {
              txt.main <- paste("Scan", scan.ids[i], "- Chromosome", chr.label)
            }
        } else {
            txt.main <- main[i]
        }
        if (!is.null(info)) {
          txt.main <- paste(txt.main, "-", info[i])
        }

        # create genotype color vector if colorGenotypes==TRUE
        cols <- .colorByGeno(colors)
        gcol <- rep(cols["NA"], length(chri))
        if (colorGenotypes) {
            genos <- getGenotype(genoData, snp=c(chr.start,chr.count), scan=c(sid,1))
            gcol[genos %in% 0] <- cols["BB"]
            gcol[genos %in% 1] <- cols["AB"]
            gcol[genos %in% 2] <- cols["AA"]
            txt.leg <- paste0(cols["AAname"], "=AA, ",
                              cols["ABname"], "=AB, ",
                              cols["BBname"], "=BB, ",
                              cols["NAname"], "=missing")
        } else {
          txt.leg <- ""
        }

        # combine plot titles for single plot
        if (type %in% c("BAF", "LRR", "R", "Theta")) {
          tmp <- paste(txt.main, txt.leg, sep="\n")
          txt.main <- tmp
          txt.leg <- tmp
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

        par(mar=c(5,4,4,2)+0.1, mgp=c(2.5,0.75,0))
        if (type == "LRR" | type == "BAF/LRR") {
            plot((posi/1e+06)[toPlot], logrratio[toPlot], xlab = "position (Mb)",
                ylab = "LRR", sub = "horizontal line = mean LRR",
                main = txt.main, ylim = c(-2, 2),
                type = "n", ...)
            for (d in 1:length(vals)) {
                abline(v = vals[d]/1e+06, col = "royalblue",
                  lty = 3, lwd = 2)
            }
            if (abst != -1)
                abline(v = abst, col = "red", lty = 2, lwd = 1.2)
            if (aben != -1)
                abline(v = aben, col = "red", lty = 2, lwd = 1.2)
            points((posi/1e+06)[toPlot], logrratio[toPlot], col=ifelse(colorGenotypes, "gray", "black"), cex=cex)
            abline(h = mninten, col = "red")
        }
        if (type == "BAF" | type == "BAF/LRR") {
            plot((posi/1e+06)[toPlot], bafreq[toPlot], type = "n",
                xlab = "position (Mb)", ylab = "BAF",
                sub = subnm, main = txt.leg, ...)
            if (abst != -1)
                abline(v = abst, col = "red", lty = 2, lwd = 1.3)
            if (aben != -1)
                abline(v = aben, col = "red", lty = 2, lwd = 1.3)
            for (d in 1:length(vals)) {
                abline(v = vals[d]/1e+06, col = "royalblue",
                  lty = 3, lwd = 2)
            }
            points((posi/1e+06)[toPlot], bafreq[toPlot], col = gcol[toPlot], cex=cex)
            abline(h = horizln, col = "gray")
        }
        if (type == "R" | type == "R/Theta") {
            plot((posi/1e+06)[toPlot], r[toPlot], xlab = "position (Mb)",
                ylab = "R", main = txt.main, type = "n", ...)
            for (d in 1:length(vals)) {
                abline(v = vals[d]/1e+06, col = "royalblue",
                  lty = 3, lwd = 2)
            }
            if (abst != -1)
                abline(v = abst, col = "red", lty = 2, lwd = 1.3)
            if (aben != -1)
                abline(v = aben, col = "red", lty = 2, lwd = 1.3)
            points((posi/1e+06)[toPlot], r[toPlot], col = gcol[toPlot], cex=cex)
        }
        if (type == "Theta" | type == "R/Theta") {
            plot((posi/1e+06)[toPlot], theta[toPlot], xlab = "position (Mb)",
                ylab = "Theta", main = txt.leg, type = "n", ...)
            for (d in 1:length(vals)) {
                abline(v = vals[d]/1e+06, col = "royalblue",
                  lty = 3, lwd = 2)
            }
            if (abst != -1)
                abline(v = abst, col = "red", lty = 2, lwd = 1.3)
            if (aben != -1)
                abline(v = aben, col = "red", lty = 2, lwd = 1.3)
            points((posi/1e+06)[toPlot], r[toPlot], col = gcol[toPlot], cex=cex)
        }

        if (ideogram) {
            par(mar=c(1,4,1,2)+0.1)
            if (ideo.zoom) {
                ideo.x <- c(min(posi[toPlot]), max(posi[toPlot]))
            } else {
                ideo.x <- c(0, lengthChromosome(chrom.char[i], "bases"))
            }
            plot(ideo.x, c(-2,2),
                 type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
            paintCytobands(chrom.char[i], units="bases", width=1, cex.leg=cex.leg)
            if (ideo.rect) rect(min(posi[toPlot]), -1.2, max(posi[toPlot]), 0.2, border="red", lwd=2)
        }
    }
}
