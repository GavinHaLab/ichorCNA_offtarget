#' plotIchorCNAzoom.R
#' author: Gavin Ha
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date:  October 8, 2019
#' description: Generate plots of copy number from On- and Off-target ichorCNA results


library(optparse)
option_list <- list(
  make_option(c("--id"), type = "character", help = "Sample ID"),
  make_option(c("--plot_funcs"), type = "character", help = "Path to file containing plotting R functions to source."),
  make_option(c("--offTargetCNFile"), type="character", help = "Path to Off-target cna.seg output file."),
  make_option(c("--onTargetCNFile"), type="character", help = "Path to On-target cna.seg output file."),
  make_option(c("--OffTargetParamFile"), type="character", help = "Path to Off-target params.txt output file."),
  make_option(c("--geneFile"), type="character", default=NULL, help = "Path to file containing list of genes with chr, start, end coordinates."),
  make_option(c("--genomeBuild"), type="character", default="hg19", help = "Genome build: hg19 or hg38. Default [%default]"),
  make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
  make_option(c("--zoom"), type = "logical", default = FALSE, help = "Zoom plot; if TRUE, then requires --chrs --start --end to be set. [Default: %default]"),
  make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to plot; string [Default: %default"),
  make_option(c("--start"), type = "integer", default = NULL, help = "Start coordinate for zoom plots"),
  make_option(c("--end"), type = "integer", default = NULL, help = "End coordinate for zoom plots"),
  make_option(c("--plotYlim"), type = "character", default = "c(-2,2)", help = "Y limits for plotting log ratio. [Default: %default]."),
  make_option(c("--plotSize"), type = "character", default = "c(5,3)", help = "width and height in inches. [Default: %default]."),
 # make_option(c("--plotFormat"), type = "character", default = "png", help = "File format of plot. E.g. pdf or png. [Default: %default]."),
  make_option(c("--outPlotFile"), type="character", help="Path to output figure file.")
 # make_option(c("--outDir"), type="character", help="Path to output directory.")
)
parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

library(data.table)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(reshape2)
library(diagram)
library(igraph)
library(tools)
library(SNPchip)
library(foreach)
library(VariantAnnotation)
library(doMC)
library(quantsmooth)

source(opt$plot_funcs)

args <- commandArgs(TRUE)
options(stringsAsFactors=F, scipen=999, bitmapType = "cairo", width=175, useDingbats = FALSE)

id <- opt$id
offTFile <- opt$offTargetCNFile
onTFile <- opt$onTargetCNFile
paramFile <- opt$OffTargetParamFile
chrStr <- as.character(eval(parse(text = "c(opt$chrs)")))
startPos <- opt$start
endPos <- opt$end
zoom <- opt$zoom
ylim <- eval(parse(text = opt$plotYlim))
geneList <- opt$geneFile
outPlotFile <- opt$outPlotFile
plotFormat <- tools::file_ext(outPlotFile)
#outDir <- opt$outDir
build <- opt$genomeBuild
genomeStyle <- opt$genomeStyle

plotSize <- eval(parse(text=opt$plotSize))
width <- plotSize[1]  #6 8
height <- plotSize[2]  #3 3.5 #4
spacing <- 3
yaxis <- "integer"
loess.span <- 0.1
numCores <- 6
coverageThres <- 0.1
minPurity <- 0.10
diffLogRThres <- 0.5
madThres <- 0.2
numPlotCols <- 3
ylimMin <- ylim[1]
ylimSV <- ylim
plotSegs <- FALSE
buffer <- 1e4
offset.factor <- 1.15 # sv drawn outside of plot
if (zoom){
  startTitle <- paste0(format(round(startPos/1e6,2), nsmall=2))
  endTitle <- paste0(format(round(endPos/1e6,2),nsmall=2), "Mb")
  xlim <- c(startPos, endPos)
  cex <- 1.5
  cytoBand <- F
  xaxt <- "s"
  plotAtCentre <- FALSE
  cnColor <- FALSE
  plotIdio <- FALSE
}else{
  xlim <- NULL
  cex <- 0.25
  cytoBand <- T
  xaxt <- "n"
  #yaxis <- "logratio"
  plotAtCentre <- FALSE
  cnColor <- FALSE
  plotIdio <- TRUE
  plotSegs <- FALSE
}
if (chrStr == "0" || is.null(chrStr) || chrStr == "None"){
  chrStr <- as.character(c(1:22, "X"))
}
seqlevelsStyle(chrStr) <- genomeStyle

if (!cnColor){
  cnCol <- rep("#000000", 30)
}else{
  cnCol <- NULL
}
cnColOff <- rep("#0072B2", 30)
cnPchOff <- 19
cnColOn <- rep("#D55E00", 30)
cnPchOn <- 17

if (!is.null(geneList) && geneList != "None"){
  genes <- read.delim(geneList, header=F, as.is=T)
}else{
  genes <- NULL
}
#midCoord <- (xlim[2] - xlim[1]) / 2 + xlim[1]
#geneAnnot <- data.frame(Gene=paste0(chrStr,":", midCoord),
#												Chr=chrStr, Start=midCoord, Stop=midCoord)
# keep track of all breakpoints #



#if (length(chrStr) > 1){
#    outPlotDir <- paste0(outDir, "/", id)
#}else{
#  outPlotDir <- outDir
#}
#dir.create(outPlotDir)
#outImage <- paste0(outDir, "/", id, ".RData")
#save.image(file=outImage)

message("Analyzing ", id)
ulpOff <- fread(offTFile)
ulpOn <- fread(onTFile)
colnames(ulpOff) <- c("Chr", "Start", "End", "copy.number", "event", "LogRatio", "subclone.status", "Corrected_Copy_Number", "Corrected_Call", "logR_Copy_Number")
colnames(ulpOn) <- c("Chr", "Start", "End", "copy.number", "event", "LogRatio", "subclone.status", "Corrected_Copy_Number", "Corrected_Call", "logR_Copy_Number")
#ulp[, state := state + 1]
#segs <- read.delim(gsub("cna.seg", "seg.txt", ulpOff), header=T, as.is=T)
#segs$median <- segs$seg.mean
#segs$chr <- segs$chrom
#segs$state <- segs$state.num
params <- read.delim(paramFile, header=T, as.is=T)
purity <- as.numeric(params[1, 2])
ploidyT <- as.numeric(params[1, 3])
normCN <- 2
ploidyS <- purity * ploidyT + (1-purity) * normCN

####### this whole block changed by Anna on 11/1/2019 ##########

if (yaxis == "integer"){
	#ulp[Chr!="X", LogRatio := log2(logRbasedCN(LogRatio, purity, ploidyT, cn=2))]
	#ulp[Chr=="X", LogRatio := log2(logRbasedCN(LogRatio, purity, ploidyT, cn=1))]

	ulpOff[Chr!="X", logR_based_copy_number := logRbasedCN(LogRatio, purity, ploidyT, cn=2)]
	ulpOff[Chr=="X", logR_based_copy_number := logRbasedCN(LogRatio, purity, ploidyT, cn=1)]

    ulpOn[Chr!="X", logR_based_copy_number := logRbasedCN(LogRatio, purity, ploidyT, cn=2)]
    ulpOn[Chr=="X", logR_based_copy_number := logRbasedCN(LogRatio, purity, ploidyT, cn=1)]


	#colName <- "logR_Copy_Number"
    colName <- "logR_based_copy_number"
}else{
	ulpOff[, LogRatio := LogRatio + log2(ploidyS / 2)]
    ulpOn[, LogRatio := LogRatio + log2(ploidyS / 2)]
	#segs$LogRatio <- segs$Median_logR
	#segs$LogRatio <- segs$LogRatio + log2(ploidyS / 2)
	colName <- "LogRatio"
}

#####################################
########## PLOT CHR RESULTS #########
#####################################
for (j in 1:length(chrStr)){
  ###################################
  if (genomeStyle == "NCBI"){
    chrTitle <- paste0("chr", chrStr[j])
  }else{
    chrTitle <- chrStr[j]
  }
  plotTitle <- paste0(id, " (", chrTitle,")")
  if (zoom){
    ylimMax <- ulpOff[Chr==chrStr[j] & Start >= xlim[1] & Start <= xlim[2], max(logR_Copy_Number, na.rm=T)] + 1
    #outPlot <- paste0(outPlotDir, "/", id, "_CNA-SV-BX_",plotType,"_chr",chrStr[j],"-",startPos,"-",endPos,".pdf")
    plotTitle <- paste0(id, " (",chrTitle,":",startTitle,"-",endTitle, ")")
  }else{
    xlim <- c(1, seqlengths(seqinfo)[chrStr[j]])
    ylimMax <- ulpOff[, max(LogRatio, na.rm=T)] + 1
  }

  ylim[2] <- min(max(ylim[2], ceiling(ylimMax)), 10)

  if (plotFormat == "png"){
    png(outPlotFile, width = width*100, height=height*100)
  }else{
    pdf(outPlotFile, width = width, height=height)
  }

  if (plotSegs) { segsToPlot <- segs } else { segsToPlot <- NULL}
  message("Plotting read depth CN")
  plotTitanIchorCNA(as.data.frame(ulpOff), segs=segsToPlot, chr=chrStr[j], colName=colName,
      cytoBand=FALSE, geneAnnot=genes, purity = purity, ploidyT = NULL, yaxis=yaxis, cnCol = cnColOff,
      yrange=ylim, xlim=xlim, spacing=spacing, xaxt=xaxt, cex = cex, gene.cex = 1, pch = cnPchOff,
      plot.title = plotTitle)

  plotTitanIchorCNA(as.data.frame(ulpOn), segs=segsToPlot, chr=chrStr[j], colName=colName,
    cytoBand=FALSE, geneAnnot=genes, purity = purity, ploidyT = NULL, yaxis=yaxis, cnCol = cnColOn,
    yrange=ylim, xlim=xlim, spacing=spacing, xaxt=xaxt, cex = cex, gene.cex = 1, pch = cnPchOn,
    plot.title = plotTitle, plot.new = FALSE)

  mtext(text = paste0("Tumor Fraction: ", format(purity, digits=3)),
            side=1, line=-1, at = xlim[2], padj=0, adj=1, cex=0.75)

  dev.off()

}

#save.image(file=outImage)
