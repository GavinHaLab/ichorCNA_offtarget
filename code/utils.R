###################################################
#################### FUNCTIONS ####################
###################################################

loadWig2GRandDT <- function(file, chrs = c(1:22, "X")){
	id <- gsub(".wig", "", basename(file))
	ulp.gr <- wigToGRanges(file)
	indChr <- orderSeqlevels(as.character(seqnames(ulp.gr)), X.is.sexchrom = TRUE)
	ulp.gr <- ulp.gr[indChr]

	# filter by chromosome
	ulp.gr <- keepSeqlevels(ulp.gr, chrs, pruning.mode="coarse")
	ulp.gr$Sample <- id

	# annotate in data.table
	ulp <- as.data.table(as.data.frame(ulp.gr))
	colnames(ulp)[c(1,6)] <- c("chr", "reads")
	ulp <- cbind(Sample = id, ulp)

	return(list(gr = ulp.gr, dt = ulp))
}

correctGCbias <- function(counts, chrs = c(1:22, "X", "Y"), gc = NULL, map = NULL, centromere = NULL,
                          flankLength = 100000, targetedSequences = NULL, genomeStyle = "NCBI", applyCorrection = TRUE,
                          mapScoreThres = 0.9, chrNormalize = c(1:22, "X", "Y")){
	seqlevelsStyle(counts) <- genomeStyle
	counts <- keepSeqlevels(counts, chrs, pruning.mode="coarse")
	if (!is.null(gc)){
		seqlevelsStyle(gc) <- genomeStyle
		counts$gc <- keepSeqlevels(gc, chrs, pruning.mode="coarse")$value
	}
	if (!is.null(map)){
		seqlevelsStyle(map) <- genomeStyle
		counts$map <- keepSeqlevels(map, chrs, pruning.mode="coarse")$value
	}
	colnames(mcols(counts))[1] <- "reads"

	# remove centromeres
	if (!is.null(centromere)){
		counts <- excludeCentromere(counts, centromere, flankLength = flankLength, genomeStyle = genomeStyle)
	}
	# keep targeted sequences
	if (!is.null(targetedSequences)){
		countsExons <- filterByTargetedSequences(counts, targetedSequences)
		counts <- counts[countsExons$ix,]
	}
	if (applyCorrection){
	## correct read counts ##
      counts <- correctReadCounts.fit(counts, chrNormalize = chrNormalize)
    }
    if (!is.null(map)) {
      ## filter bins by mappability
      counts$cor <- filterByMappabilityScore(counts$cor, map=map, mapScoreThres = mapScoreThres)
    }

  return(list(counts = counts$cor, fit.final = counts$final))
}

##################################################
###### FUNCTION TO CORRECT GC/MAP BIASES ########
##################################################
correctReadCounts.fit <- function(x, chrNormalize = c(1:22), mappability = 0.9, samplesize = 50000,
    verbose = TRUE) {
  if (length(x$reads) == 0 | length(x$gc) == 0) {
    stop("Missing one of required columns: reads, gc")
  }
  chrInd <- as.character(seqnames(x)) %in% chrNormalize
  if(verbose) { message("Applying filter on data...") }
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- quantile(x$reads[x$valid & chrInd], prob = c(0, 1 - routlier), na.rm = TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid & chrInd], prob = c(doutlier, 1 - doutlier), na.rm = TRUE)
  if (length(x$map) != 0) {
    x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
      x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  } else {
    x$ideal[!x$valid | x$reads <= range[1] |
      x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  }

  if (verbose) { message("Correcting for GC bias...") }
  set <- which(x$ideal & chrInd)
  select <- sample(set, min(length(set), samplesize))
  rough = loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final.gc = loess(predict(rough, i) ~ i, span = 0.3)
  x$cor.gc <- x$reads / predict(final.gc, x$gc)

  if (length(x$map) != 0) {
    if (verbose) { message("Correcting for mappability bias...") }
    coutlier <- 0.01
    range <- quantile(x$cor.gc[which(x$valid & chrInd)], prob = c(0, 1 - coutlier), na.rm = TRUE)
    set <- which(x$cor.gc < range[2] & chrInd)
    select <- sample(set, min(length(set), samplesize))
    final = approxfun(lowess(x$map[select], x$cor.gc[select]))
    x$cor.map <- x$cor.gc / final(x$map)
  } else {
    x$cor.map <- x$cor.gc
  }
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  return(list(cor=x, final=final.gc))
}


## note this function will write directly to ulp data.table
getRegionMeanCov <- function(ulp.gr, regions, method = "median", colToUse = "reads"){
	# find bins overlapping region(s) of interest
	#ulp <- copy(ulp)
	hits <- findOverlaps(query = ulp.gr, subject = regions)
	regionQInd <- queryHits(hits)
	regionSInd <- subjectHits(hits)

	# find mean of bins within region
	ulp <- as.data.table(as.data.frame(ulp.gr))
	ulp[regionQInd, region := regions[regionSInd, ]$Gene]
	ulp.region.mean <- na.omit(ulp[, lapply(.SD, FUN=match.fun(method), na.rm=T), by=c("Sample", "region"), .SDcols=colToUse])
	ulp.region.mean <- dcast(ulp.region.mean, Sample~region, value.var=colToUse)

	return(ulp.region.mean)
}

## find on and off target bins
getOnOffTargetBinIndex <- function(ulp.gr, baits){
	# find bins in off-target regions (i.e. not overlapping baits)
	#ulp <- copy(ulp)
	hits <- findOverlaps(query = ulp.gr, subject = baits, type = "any")
	onTargetInd <- unique(queryHits(hits))
	offTargetInd <- setdiff(1:length(ulp.gr), onTargetInd)
	return(list(onTarget = onTargetInd, offTarget = offTargetInd))
}

logRbasedCN <- function(x, purity, ploidyT, cn = 2){
	ct <- (2^x * (cn * (1 - purity) + purity * ploidyT * (cn / 2)) - cn * (1 - purity)) / purity
	ct <- sapply(ct, max, 1/2^6)
	return(ct)
}

plotBias.GC <- function(correctOutput, points = 10000, ...){
  counts <- correctOutput$counts
  #fit.rough <- correctOutput$fit.rough
  fit.final <- correctOutput$fit.final
  set <- which(counts$ideal)
  select <- sample(set, min(length(set), points))
  plot(counts$gc[select], counts$reads[select],
       col = densCols(counts$gc[select], counts$reads[select]),
       ylab = "Uncorrected Readcount", xlab = "GC content",
       ...)
  i <- seq(0, 1, by = 0.001)
  #y <- predict(fit.rough, i)
  #lines(i, y, type="l", lwd=1, col="red")
  y <- predict(fit.final, i)
  lines(i, y, type="l", lwd=1, col="red")
  coutlier = 0.001
  range <- quantile(counts$cor.gc[counts$ideal],
                    prob = c(0, 1 - coutlier), na.rm = TRUE)
  valid <- which(counts$cor.gc >= range[1] & counts$cor.gc <=
                   range[2])
  select <- intersect(valid, select)
  plot(counts$gc[select], counts$cor.gc[select],
       col = densCols(counts$gc[select], counts$cor.gc[select]),
       ylab = "GC-corrected Readcount", xlab = "GC content",
       ...)
}

# x = Data.Table containing a column called "reads"
plotReadDensity <- function(on, off, col=c("blue","red"), fill=c("blue","red")){
  gp <- ggplot(on, aes(x=reads)) +
    geom_density(color=col[1], fill=fill[1], alpha=0.5) +
    geom_density(data=off, aes(x=reads), color=col[2], fill=fill[2], alpha=0.5) +
    scale_x_continuous(trans="log10") +
    ylab("Density") + xlab("Number of reads") +
    theme_bw() +
    theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title=element_text(size=20),
          plot.title=element_text(size=24, face="bold"),
          legend.position = "none")
  return(gp)
}

plotCNlogRByChr <- function(dataIn, segs, param = NULL, colName = "copy", plotSegs = TRUE, chr=NULL, ploidy = NULL, geneAnnot=NULL, yrange=c(-4,6), xlim=NULL, col="black", xaxt = "n", cex = 0.5, gene.cex = 0.5, plot.title = NULL, spacing=4, cytoBand=T, alphaVal=1, main){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  alphaSubcloneVal <- ceiling(alphaVal / 2 * 255); class(alphaVal) = "hexmode"
  cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 26))
  subcloneCol <- c("#00FF00")
  cnCol <- paste(cnCol,alphaVal,sep="")
  names(cnCol) <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  # adjust for ploidy #
  if (!is.null(ploidy)){
    dataIn[, colName] <- as.numeric(dataIn[, colName]) + log2(ploidy / 2)
  }

  if (!is.null(chr)){
    for (i in chr){

      dataByChr <- dataIn[dataIn[,"chr"]==as.character(i),]

      #plot the data
      #if (outfile!=""){ pdf(outfile,width=10,height=6) }
      par(mar=c(spacing,8,4,2))
      #par(xpd=NA)
      coord <- (as.numeric(dataByChr[,"end"]) + as.numeric(dataByChr[,"start"]))/2
      if (is.null(xlim)){
        xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"start"]))
        xaxt <- "n"
      }
      if (is.null(plot.title)){
        plot.title <- paste("Chromosome ",i,sep="")
      }
      ## plot logR for bins ##
      plot(coord,as.numeric(dataByChr[, colName]),col=col,
           pch=16, ylim=yrange,
           xlim=xlim, xaxt = xaxt, xlab="",ylab="Copy Number\n(log2 ratio)",
           cex.lab=1.5,cex.axis=1.5, cex=cex,las=1)
      title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)
      ## plot centre line ##
      lines(c(1,xlim[2]),rep(0,2),type="l",col="grey",lwd=0.75)

      if (!is.null(geneAnnot)){
        #par(xpd=F)
        colnames(geneAnnot) <- c("Gene","Chr","Start","Stop")
        geneAnnot <- geneAnnot[geneAnnot[,"Chr"]==chr,]
        for (g in 1:dim(geneAnnot)[1]){
          print(geneAnnot[g,"Gene"])
          abline(v=as.numeric(geneAnnot[g,"Start"]),col="black",lty=3,xpd=F)
          abline(v=as.numeric(geneAnnot[g,"Stop"]),col="black",lty=3,xpd=F)
          atP <- (as.numeric(geneAnnot[g,"Stop"]) - as.numeric(geneAnnot[g,"Start"]))/2 + as.numeric(geneAnnot[g,"Start"])
          if (atP < dataByChr[1,"start"]){ atP <- dataByChr[1,"start"] }
          else if (atP > dataByChr[dim(dataByChr)[1],"start"]){ atP <- dataByChr[dim(dataByChr)[1],"start"] }
          mtext(geneAnnot[g,"Gene"],side=3,line=0,at=atP,cex=gene.cex)

        }
      }
    }
  }else{  #plot for all chromosomes
    par(mar=c(spacing,8,2,2))
    #midpt <- (as.numeric(dataIn[,"end"]) + as.numeric(dataIn[,"start"]))/2
    #coord <- getGenomeWidePositions(dataIn[,"chr"],midpt)
    coord <- getGenomeWidePositions(dataIn[,"chr"],dataIn[,"end"])
    plot(coord$posns,as.numeric(dataIn[, colName]),
         col=cnCol[as.character(dataIn[,"event"])],pch=16,xaxt="n", ylim=yrange,
         xlim=c(1,as.numeric(coord$posns[length(coord$posns)])),
         xlab="",ylab="Copy Number (log2 ratio)",
         cex.lab=1.5,cex.axis=1.5,cex=0.5,las=1,bty="n",
         #main=dataIn[1,"sample"])
         main=main)
    #plot segments
    coordEnd <- getGenomeWidePositions(segs[, "chr"], segs[, "end"])
    if (plotSegs){
      coordStart <- coordEnd$posns - (segs[, "end"] - segs[, "start"] + 1)
      xlim <- as.numeric(c(1, coordEnd$posns[length(coordEnd$posns)]))
      #col <- cnCol[as.numeric(segs[, "state"] + 1)]
      col <- col
      #write.table(segs, "~/Documents/multisample/segs_debug.seg", quote=F, sep="\t", row.names=F)  ## debug
      value <- as.numeric(segs[, "median"])
      sc.status <- as.logical(segs[, "subclone.status"])
      mat <- as.data.frame(cbind(coordStart, coordEnd$posns, value, sc.status, col))
      rownames(mat) <- 1:nrow(mat)
      ## clonal CN
      ind <- mat$sc.status == FALSE
      apply(mat[ind, ], 1, function(x){
        lines(x[1:2], rep(x[3], 2), col = x[5], lwd = 3)
        invisible()
      })
      ## subclonal CN
      if (sum(!ind) > 0){
        apply(mat[!ind, ], 1, function(x){
          lines(x[1:2], rep(x[3], 2), col = subcloneCol, lwd = 3)
          invisible()
        })
      }
    }
    lines(as.numeric(c(1,coord$posns[length(coord$posns)])),rep(0,2),type="l",col="grey",lwd=2)
    plotChrLines(dataIn[,"chr"],coordEnd$chrBkpt,yrange)
  }
}


## format of dataIn is output from /home/unix/gavinha/software/code/git/scripts/titan/analysis/combineTITAN-ichor.R
plotTitanIchorCNA <- function(dataIn, param = NULL, colName = "LogRatio", callColName="Corrected_Call", segs=NULL, chr=NULL, purity = NULL, ploidyT = NULL, geneAnnot=NULL, yrange=c(-4,6), yaxis = "logRatio", xlim=NULL, xaxt = "n", cex = 0.5, gene.cex = 0.5, pch = 16, plot.title = NULL, cnCol = NULL, spacing=4, cytoBand=T, alphaVal=1, plot.new = TRUE, ...){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  alphaSubcloneVal <- ceiling(alphaVal / 2 * 255); class(alphaVal) = "hexmode"
  subcloneCol <- c("#00FF00")
  if (is.null(cnCol)){
    cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 26))
    cnCol <- c(cnCol, "HET"="#0000FF", "DLOH"="#006400", "NLOH"="#0000FF", "ALOH"="#FF0000", "ASCNA"="#FF0000", "BCNA"="#FF0000", "UBCNA"="#FF0000")
  }else{
    cnCol.col <- as.character(cnCol[1])
    cnCol <- c(cnCol, "HET"=cnCol.col, "DLOH"=cnCol.col, "NLOH"=cnCol.col, "ALOH"=cnCol.col, "ASCNA"=cnCol.col, "BCNA"=cnCol.col, "UBCNA"=cnCol.col)
  }
  names(cnCol)[1:30] <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  #cnCol <- paste(cnCol,alphaVal,sep="")
  # adjust for ploidy #
  normCN <- 2
  if (!is.null(ploidyT) & yaxis != "integer"){
    ploidyS <- purity * ploidyT + (1-purity) * normCN
    dataIn[, colName] <- as.numeric(dataIn[, colName]) + log2(ploidyS / 2)

    if (!is.null(segs)){
      segs[, colName] <- segs[, colName] + log2(ploidyS / 2)
    }
  }


  if (!is.null(chr)){
    for (i in chr){
      dataByChr <- dataIn[dataIn[,"Chr"]==as.character(i),]
       ## set y axis labels as either integer or logR copy number
      #avgTumPloidy <- round(ploidyT)

      zero <- 0.5
      cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      #ploidyToUse <- ploidyS
      if (i == "X"){
        normCN <- 1
        zero <- 0.25
        cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      }
      if (yaxis == "integer"){
        y.ticks <- log2(cn)
        y.ticks[1] <- log2(zero)
        yrange[1] <- y.ticks[1]
        ylab <- "Copy Number"
        #dataByChr[, colName] <- log2(logRbasedCN(dataByChr[, colName], purity, ploidyT, cn=normCN))
        dataByChr[, colName] <- log2(dataByChr[, colName])
        if (!is.null(segs)){
          segs[, colName] <- log2(segs[, colName])# + log2(ploidyS / 2)
        }
        centreLine <- log2(normCN)
      }else{
        #dataByChr[, colName] <- dataByChr[, colName] + log2(ploidyS / 2)
        cnLog <- log2(cn[-which(cn==3)] / normCN)
        cn <- seq(-2,yrange[2],2)#c(-2, cn)
        y.ticks <- cn
        ylab <- "Copy Number (log2 ratio)"
        centreLine <- 0
      }

      #plot the data
      #if (outfile!=""){ pdf(outfile,width=10,height=6) }
      par(mar=c(spacing,8,4,2))
      #par(xpd=NA)
      coord <- (as.numeric(dataByChr[,"End"]) + as.numeric(dataByChr[,"Start"]))/2
      if (is.null(xlim)){
        xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"Start"]))
        xaxt <- "n"
      }
      if (is.null(plot.title)){
        plot.title <- paste("Chromosome ",i,sep="")
      }
      ## plot logR for bins ##
      if (plot.new){
        plot(coord,as.numeric(dataByChr[, colName]),col=cnCol[dataByChr[,callColName]],
             pch=pch, ylim=yrange, yaxt="n",
             xlim=xlim, xaxt = xaxt, xlab="",ylab=ylab,
             cex.lab=1.5,cex.axis=1.5, cex=cex,las=1, ...)
      }else{
        points(coord,as.numeric(dataByChr[, colName]),col=cnCol[dataByChr[,callColName]],
             pch=pch, cex=cex, ...)
      }
      axis(2, at=y.ticks, labels=cn, las=2, cex.axis=1.5)
      title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)
      ## plot centre line ##
      lines(c(1,tail(na.omit(dataByChr[,3]), 1)),rep(centreLine,2),type="l",col="grey",lwd=0.75)
      if (!is.null(segs)){
        segsByChr <- segs[segs[,"Chromosome"]==as.character(i),,drop=FALSE]
        #ind <- segsByChr$subclone.status == FALSE
        apply(segsByChr, 1, function(x){
          lines(x[c("Start","End")], rep(x[colName], 2), col = cnCol[x[callColName]], lwd = 3)
          invisible()
        })
        #if (sum(!ind) > 0){
        #  apply(segsByChr[!ind, ], 1, function(x){
        #    lines(x[c("Start","End")], rep(x["Median_logR"], 2), col = subcloneCol, lwd = 3)
        #    invisible()
        #  })
        #}
      }

      if (cytoBand==TRUE){
        require(quantsmooth)
        par(xpd = NA)
        #paintCytobands(chrom=chr, units="bases", pos=c(0,(yrange[1]-0.5)), width=0.75, legend=F)
      }

      if (!is.null(geneAnnot)){
        #par(xpd=F)
        colnames(geneAnnot) <- c("Gene","Chr","Start","Stop")
        geneAnnot <- geneAnnot[geneAnnot[,"Chr"]==chr,]
        if (nrow(geneAnnot) > 0){
        for (g in 1:dim(geneAnnot)[1]){
          print(geneAnnot[g,"Gene"])
          abline(v=as.numeric(geneAnnot[g,"Start"]),col="black",lty=3,xpd=F)
          abline(v=as.numeric(geneAnnot[g,"Stop"]),col="black",lty=3,xpd=F)
          atP <- (as.numeric(geneAnnot[g,"Stop"]) - as.numeric(geneAnnot[g,"Start"]))/2 + as.numeric(geneAnnot[g,"Start"])
          if (atP < dataByChr[1,"Start"]){ atP <- dataByChr[1,"Start"] }
          else if (atP > dataByChr[dim(dataByChr)[1],"Start"]){ atP <- dataByChr[dim(dataByChr)[1],"Start"] }
          mtext(geneAnnot[g,"Gene"],side=3,line=0,at=atP,cex=gene.cex)
          }
        }
      }
    }
  }
}


plotGenomeWide <- function(on, off=NULL, colName = "reads", plotType = "panels", geomType = "bar",
	geneAnnot=NULL, seqinfo=NULL, chrToPlot=NULL, alpha = 1, cex=0.8, ylim=NULL, xlim=NULL,
	xlab="", ylab="", col="blue", plot.title=NULL){

  colors <- c("red","blue", col); names(colors) <- c(paste0(colName,".off"), paste0(colName,".on"), colName)
  colNames <- c("Sample", "chr", "start", "end", colName)
  binSize <- on[1, end - start]
  if (is.null(xlim) && plotType != "genomeWide" && plotType != "panels"){
    xlim <- c(1, seqlengths(seqinfo[chrToPlot]))
    xaxt <- "n"
  }
  if (is.null(ylim)){
    ylim <- c(on[,min(get(colName))], on[,max(get(colName))])
  }
  if (is.null(plot.title)){
    plot.title <- paste("Chromosome ",chrToPlot,sep="")
  }

  if (!is.null(off)){
    mat <- merge(on[, mget(colNames)], off[, mget(colNames)], by=colNames[-5], all=T, suffix=c(".on",".off"))
    mat <- melt.data.table(mat, id=colNames[-5])
  }else{
    mat <- copy(on)
    mat[, value := get(colName)]
    mat[, variable := colName]
  }
  mat <- mat[chr %in% chrToPlot]
  mat <- mat[!is.na(value)]

  if (plotType == "genomeWide"){
    coord <- getGenomeWidePositions(mat$chr, mat$end, seqinfo)
  }else{
    midPts <- mat[, start + (end - start) / 2]
    coord <- NULL
    coord$posns <- midPts
  }
  mat[, positions := coord$posns]


  gp <- ggplot(mat, aes(x=positions, y=value, colour=variable, fill=variable)) +
  scale_y_continuous(limits=ylim) +
  scale_fill_manual(values=colors) + scale_color_manual(values=colors) +
  xlab(xlab) + ylab(ylab) + ggtitle(plot.title) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(face="bold", size=14),
        #panel.grid.major = element_line(colour = "grey75", linetype="dashed"),
        panel.grid.minor = element_blank(), axis.title=element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "none")
  if (geomType == "bar"){
    gp <- gp + geom_col(width = binSize/4)
  }else if (geomType == "point"){
    gp <- gp + geom_point(alpha=alpha, size=cex)
  }else{
    stop("plotGenomeWide: need to specify bar or point for geomType.")
  }

  if (plotType == "genomeWide"){ # plot genome wide and panels #
    numLines <- length(coord$chrBkpt)
    mid <- (coord$chrBkpt[1:(numLines - 1)] + coord$chrBkpt[2:numLines]) / 2
    gp <- gp + scale_x_continuous(name = "Chromosome", labels = chrToPlot, breaks = mid, limits = c(coord$chrBkpt[1], tail(coord$chrBkpt ,1))) +
      geom_vline(xintercept = coord$chrBkpt, linetype = "dotted")
  }else if (plotType == "panels"){ # additional chr panel attributes #
    gp <- gp + facet_wrap( ~ chr, ncol = 4, scales = "free_x") +
      scale_x_continuous(breaks = NULL) + expand_limits(x = 0)
  }else{ # plot for only one chr #
    gp <- gp + scale_x_continuous(name = paste0("Chromosome ", chrToPlot), breaks=NULL, limits = xlim)
  }

  return(gp)
}

ggMultiPlot <- function(..., plotlist=NULL, ncol=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)
  plotNames <- names(plots)
  # If layout is NULL, then use 'ncol' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of ncol
    layout <- matrix(seq(1, ncol * ceiling(numPlots/ncol)),
                     ncol = ncol, nrow = ceiling(numPlots/ncol), byrow = TRUE)
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
