#####
# Function to plot anomalies with stats
#####

#-------------------------------------------------------
# Accessory function

# input will be zoom type - zl for left breakpoint, zr for right breakpoint, zb for both
zoom2 <- function(pos, 	# position vector over all snps in the netCDF
	scdat, # anomaly data.frame
	j, 	# row of anomaly data.frame to be plotted
	xlim.all, 
	ztype, 	# zoom type zl=left, zr=right, zb=both
	z		# zoom factor or window size
	){
	mb <- 1000000
	leftx <- pos[scdat$left.index[j]]/mb
	rightx <- pos[scdat$right.index[j]]/mb
	dif <- rightx-leftx
	if(ztype=="zb") {left <- max(leftx-z*dif, xlim.all[1]); right <- min(rightx+z*dif, xlim.all[2])}
	if(ztype=="zl") {left <- max(leftx-z*dif, xlim.all[1]); right <- min(leftx+z*dif, xlim.all[2])}
	if(ztype=="zr") {left <- max(rightx-z*dif, xlim.all[1]); right <- min(rightx+z*dif, xlim.all[2])}			
	xlim <- c(left,right)
	return(xlim)
}	
#-------------------------------------------------------
# Main function

anomStatsPlot <- function(intenData, genoData,
	anom.stats, # data.frame with names c("anom.id","scanID", "chromosome", "left.index","right.index","method","nmark.all","nmark.elig","left.base","right.base",
		#	"nbase","non.anom.baf.med","non.anom.lrr.med","anom.baf.dev.med","anom.baf.dev.5","anom.lrr.med","nmark.baf","nmark.lrr") 
		#     this is produced by anom.seg.stats function
	snp.ineligible,	# vector of "int.id"s for intensity-only, failed snps, XTR and HLA regions
	plot.ineligible=FALSE,  # whether or not to include ineligible points in the plot
	centromere=NULL,	# if a data.frame of centromere positions is given here, they will be shown as a gray rectangle on the plot
				# centromere is a data.frame with variables "chromosome","left.base","right.base"
	brackets = c("none","bases","markers"),	# re plotting brackets around breakpoints - none, use base length, use number of markers (note that using markers give asymmetric brackets)
	brkpt.pct = 10,  # percent of anomaly length in bases (or number of markers) for width of brackets
	whole.chrom=FALSE,  # logical to plot the whole chromosome or not (over-rides win and zoom arguments)
	win=5,			# size of the window surrounding the anomaly to plot = a multiple of anomaly length
	win.calc=FALSE,		# calculate window size from anomaly length; over-rides win and gives window of fixed length given by win.fixed
	win.fixed = 1,		# number of megabases for window size when win.calc=TRUE
	zoom=c("both","left","right"), # indicates whether plot includes the whole anomaly ("both") or zooms on just the left or right breakpoint; "both" is default
	info=NULL,		# character vector of extra information to include in the main title of the upper plot                          
	main.txt = NULL,                          
	type = c("BAF/LRR", "BAF", "LRR"),
	cex=0.5		# cex value for points on the plots                       
){

	# v3 adds option to plot vertical lines to bracket each breakpoint, width of bracket being a percentage of the length or number of markers 
	# 	in the anomaly
	#     also uses new var names for 'anom' data.frame (from the new anom.seg.stats function)
	#	allow for missing value codes other than -1 for genotypes
	# v4 adds options to identify the centromere as a shaded rectangle and to make whole-chromosome plot
	# v5 adds option to add a second row of identifier information to the main title
	# v6 adds arguments 'zoom','win.calc','win.fixed'; fix use of to bl.ncdf.file and geno.ncdf.file; add plotting of terminal markers as orange vertical lines (among all markers, intensity-only or not)
	# v7 fix arg order and defaults; when plot.ineligible=T these will be plotted only on LRR, not on BAF
	# v8 add anom.id to plot title
	# v9 when plotting ineligible points for lrr, put them underneath the eligible


  # check that intenData has BAF/LRR
  if (!hasBAlleleFreq(intenData)) stop("BAlleleFreq not found in intenData")
  if (!hasLogRRatio(intenData)) stop("LogRRatio not found in intenData")
  
  # check that dimensions of intenData and genoData are equal
  intenSnpID <- getSnpID(intenData)
  genoSnpID <- getSnpID(genoData)
  if (!all(intenSnpID == genoSnpID)) stop("snp dimensions of intenData and genoData differ")
  intenScanID <- getScanID(intenData)
  genoScanID <- getScanID(genoData)
  if (!all(intenScanID == genoScanID)) stop("scan dimensions of intenData and genoData differ")
  
	intid <- intenSnpID
	sid <- intenScanID
	indices <- 1:length(intid)
	chrom <- getChromosome(intenData)
	pos <- getPosition(intenData)

	# check centromere
	if(!is.null(centromere)){
		if(!class(centromere)=="data.frame") stop("centromere must be a data.frame")
		if(!all(is.element(c("chrom","left.base","right.base"), names(centromere)))) stop("names of centromere must include chrom, left.base, right.base")
                centromere$chrom[is.element(centromere$chrom, "X")] <- XchromCode(intenData)
                centromere$chrom[is.element(centromere$chrom, "Y")] <- YchromCode(intenData)
                centromere$chrom[is.element(centromere$chrom, "XY")] <- XYchromCode(intenData)
                centromere$chrom <- as.integer(centromere$chrom)
		centromere <- centromere[, c("chrom","left.base","right.base")]
		if(!all(unlist(lapply(centromere, class))=="integer")) stop("required centromere variables must be of class integer")
		if(!all(is.element(c(1:22, XchromCode(intenData)),centromere$chrom))) stop("centromere positions must be given for chroms 1-22 and X")
	}
	cent.col <- "lightblue1"
	term.col <- "orange"

	# zoom
	zoom <- match.arg(zoom)

	# check anom.stats
	if(class(anom.stats)!="data.frame") stop("anom.stats should be a data.frame")
	anames <- c("anom.id","scanID", "chromosome", "left.index","right.index","method","nmark.all","nmark.elig","left.base","right.base",
		"nbase","non.anom.baf.med","non.anom.lrr.med","anom.baf.dev.med","anom.baf.dev.5","anom.lrr.med","nmark.baf","nmark.lrr")
	if(!all(is.element(anames,names(anom.stats)))) stop("anom.stats does not have required variable names")
	if(!all(is.element(c(anom.stats$left.index,anom.stats$right.index),indices))) stop("left and/or right are not within range of snp indices")

	# check snp.ineligible
	if(!all(is.element(snp.ineligible,intid))) stop("snp.ineligible are not int.ids")

	# check win
	if(!is.null(win)) if(!(is.numeric(win) & length(win==1) & win>0)) stop("win must be null or a positive number")
	if(win.calc==FALSE){
		if(is.null(win)) stop("specify win when win.calc=F")
		} 
	
	# check brackets and brkpt.pct
        brackets <- match.arg(brackets)
	if(!(brkpt.pct>0 & brkpt.pct<50)) stop("brkpt.pct must be >0 and <50")

	# check info
	if(!is.null(info)) {
		if(!(length(info)==nrow(anom.stats) & is.character(info))) stop("info should be a character vector of length equal to rows of anom.stats")
	}

	# logical for ineligible
	ne <- is.element(intid, snp.ineligible)
	
	# bases per megabase
	mb <- 1000000

        type <- match.arg(type)
        if (type == "BAF/LRR") {
          par(mfrow=c(2,1))
        }
	for(i in 1:nrow(anom.stats)) {
		sind <- which(is.element(sid, anom.stats$scanID[i]))
		chr <- chrom==anom.stats$chromosome[i]
		left <- anom.stats$left.index[i]
		right <- anom.stats$right.index[i]
		baf <- getBAlleleFreq(intenData, snp=c(1,-1), scan=c(sind,1))
		lrr <- getLogRRatio(intenData, snp=c(1,-1), scan=c(sind,1))
		geno <- getGenotype(genoData, snp=c(1,-1), scan=c(sind,1))

		#check base positions
		chk <- pos[left]==anom.stats$left.base[i] & pos[right]==anom.stats$right.base[i]
		if(!all(chk)) stop(paste("left.base and/or right.base in anom.stats table row", i, "not correct"))

		# calculate window size if indicated
		if(win.calc){
			alen <- anom.stats[i,"nbase"]/mb	# anomaly length in Mb
			win <- win.fixed/alen		
		}	

		# make missing genotype calls equal -1 (regardless of missing value code)
		geno[is.na(geno)] <- -1

		# color-coding for genotype calls
		gcol <- rep(NA, length(geno))
		gcol[geno==2] <- "red"
		gcol[geno==1] <- "green"
		gcol[geno==0] <- "blue"
		gcol[geno== -1] <- "black"
		gcol[ne] <- "pink"
		# color-coding for lrr
		lcol <- rep("gray", length(lrr))
		lcol[ne] <- "pink"
		
		# set xlim values for base position, leftp and rightp
		# terminal positions for the selected chromosome
			posc <- pos[chr]
			n <- length(posc)
			xlim.all <- c(posc[1], posc[n])/mb
		# anom length
		if(whole.chrom==TRUE) { xlim <- xlim.all	# set xlim to be the start and end of the chromosome
		} else {
			if(zoom=="left")  xlim <- GWASTools:::zoom2(pos=pos, scdat=anom.stats, j=i, xlim.all=xlim.all, ztype="zl", z=win)	# zoom2 accessory function defined above
			if(zoom=="right") xlim <- GWASTools:::zoom2(pos=pos, scdat=anom.stats, j=i, xlim.all=xlim.all, ztype="zr", z=win)
			if(zoom=="both")  xlim <- GWASTools:::zoom2(pos=pos, scdat=anom.stats, j=i, xlim.all=xlim.all, ztype="zb", z=win)
		}				


		# set lower limit on lrr to be -2 or the min within the anomaly, whichever is smaller
		lrr.min <- min(lrr[left:right], na.rm=T)
		lrr.min <- min(lrr.min,-2)

		# select points to plot
		if(plot.ineligible==FALSE) {sel.lrr <- chr & !ne
			}  else {sel.lrr <- chr}
		sel.baf <- chr & !ne

		# brackets around breakpoints
		if(brackets=="bases"){
			br <- (brkpt.pct/100)*(anom.stats$nbase[i])
			left.lower <- pos[left]-br
			left.upper <- pos[left]+br
			right.lower <- pos[right]-br
			right.upper <- pos[right]+br
			bkt.pos <- c(left.lower/mb, left.upper/mb, right.lower/mb, right.upper/mb)
		}
		if(brackets=="markers"){
			br <- ceiling((brkpt.pct/100)*(anom.stats$nmark.elig[i])) # number of eligible markers to include in bracket
			elig <- !ne
			bre <- 0; j <- 0; while(bre<br) { j <- j+1; bre <- bre + elig[left-j]}
			left.lower <- pos[left-j]
			table(elig[(left-j):(left-1)])
			bre <- 0; j <- 0; while(bre<br) { j <- j+1; bre <- bre + elig[left+j-1]}  # includes the breakpoint marker
			left.upper <- pos[left+j-1]
			table(elig[left:(left+j-1)])
			bre <- 0; j <- 0; while(bre<br) { j <- j+1; bre <- bre + elig[right-j]}
			right.lower <- pos[right-j]
			table(elig[(right-j):(right-1)])
			bre <- 0; j <- 0; while(bre<br) { j <- j+1; bre <- bre + elig[right+j-1]} # includes the breakpoint marker
			right.upper <- pos[right+j-1]
			table(elig[right:(right+j-1)])
			bkt.pos <- c(left.lower/mb, left.upper/mb, right.lower/mb, right.upper/mb)
		}

		# get centromere start and end points
		if(!is.null(centromere)){
			cent <- centromere[centromere$chrom==anom.stats$chromosome[i],][1,] # if >1 centromere position given, the first is used
			if(nrow(cent)>1) stop(paste("more than one row in centromere for chromosome",i))
		}
			
		# plot LRR
        if (type == "LRR" | type == "BAF/LRR") {
          if (is.null(main.txt)) {
		mtxt <- paste("anom",anom.stats$anom.id[i],"-  snum", anom.stats$scanID[i], "- chrom", anom.stats$chrom[i], "-", anom.stats$sex[i], "-", anom.stats$method[i])
              } else {
                mtxt <- main.txt
              }
		if(brackets!="none") mtxt <- paste(mtxt, "- gray brackets =", brkpt.pct, "% of", brackets)
		if(!is.null(info)) mtxt <- paste(mtxt, "\n", info[i])
		stxt <- "red=AA, green=AB, blue=BB, pink=ineligible, black=other missing"
		plot(pos[sel.lrr]/mb, lrr[sel.lrr], xlab="position (Mb)", ylab="LRR", type="n", ylim=c(lrr.min,2), xlim=xlim, main=mtxt) 
		if(!is.null(centromere)){
			rect(cent$left.base/mb,lrr.min,cent$right.base/mb,2, col=cent.col, border=cent.col)
		}
		abline(v=xlim.all, col=term.col)  # terminal markers
		abline(h=0, col="gray")
		abline(v=c(pos[left]/mb,pos[right]/mb), col="black")
		if(brackets!="none") abline(v=bkt.pos, col="gray")
		if(plot.ineligible) points(pos[chr & ne]/mb, lrr[chr & ne], col="pink", cex=cex)
		points(pos[chr & !ne]/mb, lrr[chr & !ne], col="darkgray", cex=cex)
		# Add anom stats
		abline(h=anom.stats$non.anom.lrr.med[i], col="red")
		tmp <- anom.stats$anom.lrr.med[i]
		segments(pos[left]/mb, tmp, pos[right]/mb, tmp, col="red", lty=2)
              }

		# plot BAF
        if (type == "BAF" | type == "BAF/LRR") {
		mtxt <- paste("red=AA, green=AB, blue=BB, pink=ineligible, black=other missing",
                              "horiz solid red = non-anom median, horiz dashed red = anom median",
                              sep="\n")
		plot(pos[sel.baf]/mb, baf[sel.baf], xlab="position (Mb)", ylab="BAF", type="n", ylim=c(0,1), xlim=xlim, main=mtxt) 
		if(!is.null(centromere)){
			rect(cent$left.base/mb,0,cent$right.base/mb,1, col=cent.col, border=cent.col)
		}
		abline(v=xlim.all, col=term.col)    # terminal markers
		abline(v=c(pos[left]/mb,pos[right]/mb), col="black")
		abline(h=c(0.5,1/3,2/3,0,1), col="pink")
		if(brackets!="none") abline(v=bkt.pos, col="gray")
		points(pos[sel.baf]/mb, baf[sel.baf], col=gcol[sel.baf], cex=cex)
		# Add anom stats
		abline(h=anom.stats$non.anom.baf.med[i], col="red")
		tmp <- anom.stats$non.anom.baf.med[i] + anom.stats$anom.baf.dev.med[i] 
		segments(pos[left]/mb, tmp, pos[right]/mb, tmp, col="red", lty=2)
		tmp <- anom.stats$non.anom.baf.med[i] - anom.stats$anom.baf.dev.med[i] 
		segments(pos[left]/mb, tmp, pos[right]/mb, tmp, col="red", lty=2)
              }
	
	}  # end loop over anomalies

} # end function definition


	
