# file:   ichorCNA_offTarget.R
# author: Gavin Ha, Ph.D.
#         Fred Hutch
# contact: <gha@fredhutch.org>

# ichorCNA: https://github.com/broadinstitute/ichorCNA
# date:   July 24, 2019
# description:

library(optparse)

option_list <- list(
  make_option(c("--TUMWIG"), type = "character", help = "Path to tumor WIG file. Required."),
  make_option(c("--NORMWIG"), type = "character", default=NULL, help = "Path to normal WIG file. Default: [%default]"),
  make_option(c("--baitBedTum"), type = "character", help = "Path to tumor bait bed file. Required"),
  make_option(c("--baitBedNorm"), type = "character", default=NULL, help = "Path to normal bait bed file. Optional. Will use --baitBedTum if not provided."),
  make_option(c("--gcWig"), type = "character", help = "Path to GC-content WIG file; Required"),
  make_option(c("--mapWig"), type = "character", default=NULL, help = "Path to mappability score WIG file. Default: [%default]"),
  make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
  make_option(c("--centromere"), type="character", default=NULL, help = "File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package. Default: [%default]"),
  make_option(c("--genomeBuild"), type="character", default="hg19", help="Geome build. Default: [%default]"),
  make_option(c("--genomeStyle"), type="character", default="NCBI", help="Genome style {NCBI, UCSC}. Default [%default]"),
  make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases. Default: [%default]"),
  make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze. Default: [%default]"),
  make_option(c("--normalizeMaleX"), type="logical", default=TRUE, help = "If male, then normalize chrX by median. Default: [%default]"),
  make_option(c("--fracReadsInChrYForMale"), type="numeric", default=1e-05, help = "Threshold for fraction of reads in chrY to assign as male. Default: [%default]"),
  	make_option(c("--plotFileType"), type="character", default="pdf", help = "File format for output plots. Default: [%default]"),
	make_option(c("--plotYLim"), type="character", default="c(-2,2)", help = "ylim to use for chromosome plots. Default: [%default]"),
  make_option(c("--outDir"), type="character", default="./", help = "Output Directory. Default: [%default]"),
  make_option(c("--libdir"), type = "character", default=NULL, help = "Script library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=5, stringsAsFactors=F)

library(HMMcopy)
library(data.table)
library(GenomicRanges)
library(ggplot2)
#library(ggbio)
library(SNPchip)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

id <- opt$id
tumour_file <- opt$TUMWIG
normal_file <- opt$NORMWIG
gcWig <- opt$gcWig
mapWig <- opt$mapWig
baitFile.tum <- opt$baitBedTum  # "0" if none specified
baitFile.norm <- opt$baitBedNorm  # "0" if none specified
centromere <- opt$centromere
flankLength <- opt$rmCentromereFlankLength
normal <- eval(parse(text = opt$normal))
normalizeMaleX <- as.logical(opt$normalizeMaleX)
fracReadsInChrYForMale <- opt$fracReadsInChrYForMale
outDir <- opt$outDir
libdir <- opt$libdir
plotFileType <- opt$plotFileType
plotYLim <- eval(parse(text=opt$plotYLim))
gender <- NULL
outImage <- paste0(outDir,"/", id,".RData")
chrs <- eval(parse(text = opt$chrs))
chrs.all <- c(chrs, "Y")
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize)))
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
seqlevelsStyle(chrs) <- genomeStyle
seqlevelsStyle(chrs.all) <- genomeStyle
seqlevelsStyle(chrNormalize) <- genomeStyle
mapScoreThres <- 0.8

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
outPlotDir <- paste0(outDir, "/", id)
dir.create(outPlotDir)

save.image(outImage)


if (!is.null(libdir)){
	source(paste0(libdir,"/R/utils.R"))
	source(paste0(libdir,"/R/segmentation.R"))
	source(paste0(libdir,"/R/EM.R"))
	source(paste0(libdir,"/R/output.R"))
	source(paste0(libdir,"/R/plotting.R"))
} else {
    library(ichorCNA)
}
#source("~/software/git/scripts/hmmcopy/analysis/offTargetWES/ichorCNA_utils.R")
#source("~/software/git/scripts/hmmcopy/analysis/offTargetWES/ichorCNA_plotting.R")

###Changed this to relative filepath --Anna 9/13/19###
source("../utils.R")

## load seqinfo
seqinfo <- getSeqInfo(genomeBuild, genomeStyle)

# load bait interval file for this tumor normal pair
# exclude patient if bait file diff for tumor and normal capture
baits <- fread(baitFile.tum, skip=0)
colnames(baits) <- c("chr", "start", "end", "name", "strands")
baits  <- as(baits , "GRanges")
seqlevelsStyle(baits) <- genomeStyle
if (!is.null(baitFile.norm)){
  baits.norm <- fread(baitFile.norm, skip=0)
  colnames(baits.norm ) <- c("chr", "start", "end", "name", "strands")
  baits.norm  <- as(baits.norm , "GRanges")
  seqlevelsStyle(baits.norm) <- genomeStyle
	if (length(baitFile.tum) != length(baitFile.norm))
		stop("Mismatch in tumor and normal bait interval files or NA in at least one.")
}

# load gc and map wigs into RangedData
gc <- wigToGRanges(gcWig)
seqlevelsStyle(gc) <- genomeStyle
gc <- keepChr(gc, chrs.all)
map <- wigToGRanges(mapWig)
seqlevelsStyle(map) <- genomeStyle
map <- keepChr(map, chrs.all)
centromere <- fread(centromere)
seqlevelsStyle(centromere$Chr) <- genomeStyle

# tumor sample
ulpT <- loadWig2GRandDT(tumour_file, chrs.all)
onOffTargetInd.tum <- getOnOffTargetBinIndex(ulpT$gr, baits)
ulpT.onTarget.gcCor <- correctGCbias(ulpT$gr[onOffTargetInd.tum$onTarget,],
	chrs = chrs, genomeStyle = genomeStyle,
	gc = gc[onOffTargetInd.tum$onTarget,], map = map[onOffTargetInd.tum$onTarget,],
	centromere = centromere, mapScoreThres = mapScoreThres, chrNormalize = chrNormalize)
ulpT.offTarget.gcCor <- correctGCbias(ulpT$gr[onOffTargetInd.tum$offTarget,],
	chrs = chrs, genomeStyle = genomeStyle,
	gc = gc[onOffTargetInd.tum$offTarget,], map = map[onOffTargetInd.tum$offTarget,],
	centromere = centromere, mapScoreThres = mapScoreThres, chrNormalize = chrNormalize)
# normal sample
ulpN <- loadWig2GRandDT(normal_file, chrs)
onOffTargetInd.norm <- getOnOffTargetBinIndex(ulpN$gr, baits)
ulpN.onTarget.gcCor <- correctGCbias(ulpN$gr[onOffTargetInd.norm$onTarget,],
	chrs = chrs, gc = gc[onOffTargetInd.norm$onTarget,], map = map[onOffTargetInd.norm$onTarget,],
	centromere = centromere, mapScoreThres = mapScoreThres, chrNormalize = chrNormalize)
ulpN.offTarget.gcCor <- correctGCbias(ulpN$gr[onOffTargetInd.norm$offTarget,],
	chrs = chrs, gc = gc[onOffTargetInd.norm$offTarget,], map = map[onOffTargetInd.norm$offTarget,],
	centromere = centromere, mapScoreThres = mapScoreThres, chrNormalize = chrNormalize)

## get gender ##
counts <- loadReadCountsFromWig(wigToGRanges(normal_file), chrs=chrs, gc=gc, map=map, fracReadsInChrYForMale =
 				fracReadsInChrYForMale, chrXMedianForMale = Inf,
 				centromere=centromere, flankLength = flankLength, targetedSequences=NULL,
				genomeStyle = genomeStyle, chrNormalize = chrNormalize, mapScoreThres = mapScoreThres)
gender <- counts$gender$gender
# if male, then just normalize chrX to median (ULP and WES)
if (gender=="male" && normalizeMaleX){
	# on target
	chrXInd.on <- grep("X", as.character(seqnames(ulpT.onTarget.gcCor$counts)))
	chrXMedian.on <- median(ulpT.onTarget.gcCor$counts$copy[chrXInd.on], na.rm = TRUE)
	ulpT.onTarget.gcCor$counts$copy[chrXInd.on] <- ulpT.onTarget.gcCor$counts$copy[chrXInd.on] - chrXMedian.on
	# off target
	chrXInd.off <- grep("X", as.character(seqnames(ulpT.offTarget.gcCor$counts)))
	chrXMedian.off <- median(ulpT.offTarget.gcCor$counts$copy[chrXInd.off], na.rm = TRUE)
	ulpT.offTarget.gcCor$counts$copy[chrXInd.off] <- ulpT.offTarget.gcCor$counts$copy[chrXInd.off] - chrXMedian.off
}

# normalize tumor by normal
ulpT.onTarget.gcCor$counts$copy.norm <- ulpT.onTarget.gcCor$counts$copy - ulpN.onTarget.gcCor$counts$copy
ulpT.offTarget.gcCor$counts$copy.norm <- ulpT.offTarget.gcCor$counts$copy - ulpN.offTarget.gcCor$counts$copy


# convert GRanges to data.table
ulpT.onTarget.gcCor.dt <- cbind(as.data.table(ulpT.onTarget.gcCor$counts))
setnames(ulpT.onTarget.gcCor.dt, c("seqnames"), c("chr"))
ulpT.offTarget.gcCor.dt <- cbind(as.data.table(ulpT.offTarget.gcCor$counts))
setnames(ulpT.offTarget.gcCor.dt, c("seqnames"), c("chr"))
ulpN.onTarget.gcCor.dt <- cbind(as.data.table(ulpN.onTarget.gcCor$counts))
setnames(ulpN.onTarget.gcCor.dt, c("seqnames"), c("chr"))
ulpN.offTarget.gcCor.dt <- cbind(as.data.table(ulpN.offTarget.gcCor$counts))
setnames(ulpN.offTarget.gcCor.dt, c("seqnames"), c("chr"))

# combine tumor and normal into same data.table
onTarget.gcCor.dt <- merge(ulpT.onTarget.gcCor.dt, ulpN.onTarget.gcCor.dt, by=c("chr","start","end","width","strand","gc","map"), suffix=c(".tumor",".normal"))
offTarget.gcCor.dt <- merge(ulpT.offTarget.gcCor.dt, ulpN.offTarget.gcCor.dt, by=c("chr","start","end","width","strand","gc","map"), suffix=c(".tumor",".normal"))
## output results to text files
outFile <- paste0(outDir, "/", id, "_onTarget_cor.txt")
fwrite(onTarget.gcCor.dt, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")
outFile <- paste0(outDir, "/", id, "_offTarget_cor.txt")
fwrite(offTarget.gcCor.dt, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

# collect metrics for reads genome-wide and chrX for tumor and normal, on- and off-target
onTarget.chrXCovT <- ulpT$dt[onOffTargetInd.tum$onTarget][grepl("X", chr), median(reads)]
offTarget.chrXCovT <- ulpT$dt[onOffTargetInd.tum$offTarget][grepl("X", chr), median(reads)]
onTarget.genomeCovT <- ulpT$dt[onOffTargetInd.tum$onTarget][, median(reads)]
offTarget.genomeCovT <- ulpT$dt[onOffTargetInd.tum$offTarget][, median(reads)]
onTarget.chrXCovN <- ulpN$dt[onOffTargetInd.norm$onTarget][grepl("X", chr), median(reads)]
offTarget.chrXCovN <- ulpN$dt[onOffTargetInd.norm$offTarget][grepl("X", chr), median(reads)]
onTarget.genomeCovN <- ulpN$dt[onOffTargetInd.norm$onTarget][, median(reads)]
offTarget.genomeCovN <- ulpN$dt[onOffTargetInd.norm$offTarget][, median(reads)]

mat <- cbind(Sample=id, Sex=gender,
		data.table(Tum.offTarget=cbind(Num.bins=length(onOffTargetInd.tum$offTarget),
							  genome.reads=offTarget.genomeCovT, chrX.reads=offTarget.chrXCovT),
		  Tum.onTarget=cbind(Num.bins=length(onOffTargetInd.tum$onTarget),
							 genome.reads=onTarget.genomeCovT, chrX.reads=onTarget.chrXCovT),
		  Norm.offTarget=cbind(Num.bins=length(onOffTargetInd.norm$offTarget),
							   genome.reads=offTarget.genomeCovN, chrX.reads=offTarget.chrXCovN),
		  Norm.onTarget=cbind(Num.bins=length(onOffTargetInd.norm$onTarget),
							  genome.reads=onTarget.genomeCovN, chrX.reads=onTarget.chrXCovN)))
outFile <- paste0(outDir, "/", id, "_readStats.txt")
fwrite(mat, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

save.image(outImage)

## plot bias
# on-target - tumor & normal
outPlot <- paste0(outPlotDir, "/", id, "_biasCor_tumor.pdf")
pdf(outPlot, width=7, height=7)
par(mfrow=c(2,2))
tum.rd.on.plot <- ulpT.onTarget.gcCor
tum.rd.on.plot$copy <- tum.rd.on.plot$copy.norm
plotBias.GC(tum.rd.on.plot, pch = 20, cex = 0.5, main="On-Target")
tum.rd.off.plot <- ulpT.offTarget.gcCor
tum.rd.off.plot$copy <- tum.rd.off.plot$copy.norm
plotBias.GC(tum.rd.off.plot, pch = 20, cex = 0.5, main="Off-Target")
dev.off()

# off-target - tumor & normal
outPlot <- paste0(outPlotDir, "/", id, "_biasCor_normal.pdf")
pdf(outPlot, width=7, height=7)
par(mfrow=c(2,2))
norm.rd.on.plot <- ulpN.onTarget.gcCor
norm.rd.on.plot$copy <- norm.rd.on.plot$copy.norm
plotBias.GC(norm.rd.on.plot, pch = 20, cex = 0.5, main="On-Target")
norm.rd.off.plot <- ulpN.offTarget.gcCor
norm.rd.off.plot$copy <- norm.rd.off.plot$copy.norm
plotBias.GC(norm.rd.off.plot, pch = 20, cex = 0.5, main="Off-Target")
dev.off()

## plot read densities
gp <- plotReadDensity(ulpT.onTarget.gcCor.dt, ulpT.offTarget.gcCor.dt)
outPlot <- paste0(outPlotDir, "/", id, "_density_tumor.pdf")
ggsave(gp, file=outPlot, width=5, height=5)
gp <- plotReadDensity(ulpN.onTarget.gcCor.dt, ulpN.offTarget.gcCor.dt)
outPlot <- paste0(outPlotDir, "/", id, "_density_normal.pdf")
ggsave(gp, file=outPlot, width=5, height=5)

## plot genome-wide raw read counts per bin ##
gp <- plotGenomeWide(ulpT.onTarget.gcCor.dt, ulpT.offTarget.gcCor.dt, chrToPlot=c(1:22,"X"), colName = "reads", plotType = "panels", geomType = "bar", geneAnnot=NULL, seqinfo=seqinfo, ylim=NULL, xlim=NULL, xlab="", ylab="Number of reads", plot.title=NULL)
outPlot <- paste0(outPlotDir, "/", id, "_readCov_tumor.pdf")
ggsave(gp, file=outPlot, width=12, height=8)
gp <- plotGenomeWide(ulpN.onTarget.gcCor.dt, ulpN.offTarget.gcCor.dt, chrToPlot=c(1:22,"X"), colName = "reads", plotType = "panels", geomType = "bar", geneAnnot=NULL, seqinfo=seqinfo, ylim=NULL, xlim=NULL, xlab="", ylab="Number of reads", plot.title=NULL)
outPlot <- paste0(outPlotDir, "/", id, "_readCov_normal.pdf")
ggsave(gp, file=outPlot, width=12, height=8)

## plot genome-wide panels GC-corrected counts and matched normal normalized
gp <- plotGenomeWide(ulpT.onTarget.gcCor.dt, chrToPlot=c(1:22,"X"), colName = "copy.norm", plotType = "panels", geomType = "point", geneAnnot=NULL, seqinfo=seqinfo, ylim=c(-4,4), xlim=NULL, xlab="", ylab="Number of reads", col="blue", plot.title=NULL)
outPlot <- paste0(outPlotDir, "/", id, "_corLogR_tumor_onTarget.pdf")
ggsave(gp, file=outPlot, width=12, height=8)
gp <- plotGenomeWide(ulpT.offTarget.gcCor.dt, chrToPlot=c(1:22,"X"), colName = "copy.norm", plotType = "panels", geomType = "point", geneAnnot=NULL, seqinfo=seqinfo, ylim=c(-4,4), xlim=NULL, xlab="", ylab="Number of reads", col="red", plot.title=NULL)
outPlot <- paste0(outPlotDir, "/", id, "_corLogR_tumor_offTarget.pdf")
ggsave(gp, file=outPlot, width=12, height=8)

## plot genome-wide (no panels, just full profile) counts and matched normal normalized
gp.on <- plotGenomeWide(ulpT.onTarget.gcCor.dt, chrToPlot=c(1:22,"X"), colName = "copy.norm", plotType = "genomeWide", geomType = "point", geneAnnot=NULL, seqinfo=seqinfo, alpha=0.5, ylim=c(-4,4), xlim=NULL, xlab="", ylab="Log2 Ratio", col="blue", plot.title="On-Target")
gp.off <- plotGenomeWide(ulpT.offTarget.gcCor.dt, chrToPlot=c(1:22,"X"), colName = "copy.norm", plotType = "genomeWide", geomType = "point", geneAnnot=NULL, seqinfo=seqinfo, alpha=0.5, ylim=c(-4,4), xlim=NULL, xlab="", ylab="Log2 Ratio", col="red", plot.title="Off-Target")
outPlot <- paste0(outPlotDir, "/", id, "_corLogR_tumor_genomeWide.pdf")
pdf(outPlot, width=10, height=5)
ggMultiPlot(plotlist=list(gp.on, gp.off), ncol=1, layout=NULL)
dev.off()


## plot per chromosome GC-corrected counts and matched normal normalized
ylim <- plotYLim
for (i in 1:length(chrs)){
	xlim <- c(1, seqlengths(seqinfo)[chrs[i]])
	outPlot <- paste0(outPlotDir, "/", id, "_corLogR_tumor_chr",chrs[i],".pdf")
	pdf(outPlot, width=10, height=6)
	chr.ideo <- chrs[i]
	seqlevelsStyle(chr.ideo) <- "UCSC"
	par(mfrow=c(2,1))
	par(xpd=NA)
	plot.title <- paste0("On-Target (", chr.ideo, ")")
	plotCNlogRByChr(as.data.frame(ulpT.onTarget.gcCor.dt), segs = NULL, plotSegs=F, colName = "copy.norm", chr=chrs[i],
						ploidy = NULL, cytoBand=F, yrange=ylim, xlim=xlim, spacing=4, col="black", plot.title=plot.title)
	pI <- plotIdiogram(chrs[i], build=genomeBuild, unit="bp", label.y=ylim[1]-abs(ylim[1])*0.8, label.cytoband=FALSE,
						 new=FALSE, ylim=c(ylim[1]-abs(ylim[1])*0.40,ylim[1]-abs(ylim[1])*0.25))
	plot.title <- paste0("Off-Target (", chr.ideo, ")")
	plotCNlogRByChr(as.data.frame(ulpT.offTarget.gcCor.dt), segs = NULL, plotSegs=F, colName = "copy.norm", chr=chrs[i],
						ploidy = NULL, cytoBand=F, yrange=ylim, xlim=xlim, spacing=4, col="black", plot.title=plot.title)
	pI <- plotIdiogram(chrs[i], build=genomeBuild, unit="bp", label.y=ylim[1]-abs(ylim[1])*0.8, label.cytoband=FALSE,
						 cex=0.1, new=FALSE, ylim=c(ylim[1]-abs(ylim[1])*0.40,ylim[1]-abs(ylim[1])*0.25))
	dev.off()

 #gp.on <- plotGenomeWide(ulpT.onTarget.gcCor.dt, chrToPlot=chrs[i], colName = "copy.norm", plotType = "chr", geomType = "point", geneAnnot=NULL, seqinfo=seqinfo, ylim=ylim xlim=NULL, xlab="", ylab="Log2 Ratio", col="blue", plot.title=NULL)
  #gp.off <- plotGenomeWide(ulpT.offTarget.gcCor.dt, chrToPlot=chrs[i], colName = "copy.norm", plotType = "chr", geomType = "point", geneAnnot=NULL, seqinfo=seqinfo, ylim=ylim, xlim=NULL, xlab="", ylab="Log2 Ratio", col="red", plot.title=NULL)

  #chr.ideo <- chrs[i]
  #seqlevelsStyle(chr.ideo) <- "UCSC"
  #gp.ideo <- plotIdeogram(subchr = chr.ideo, genome = genomeBuild)
  #plot.new()
  #ggMultiPlot(plotlist=list(gp.on, gp.off), ncol=1, layout=NULL)

# chr.ideo <- chrs[i]
# seqlevelsStyle(chr.ideo) <- "UCSC"
# data1.x <- ulpT.onTarget.gcCor.dt[chr==chrs[i], start]
# data1.y <- ulpT.onTarget.gcCor.dt[chr==chrs[i], copy.norm]
# data2.x <- ulpT.offTarget.gcCor.dt[chr==chrs[i], start]
# data2.y <- ulpT.offTarget.gcCor.dt[chr==chrs[i], copy.norm]
# pp <- getDefaultPlotParams(plot.type=1)
# pp$topmargin <- 25
# pp$bottommargin <- 25
# pp$ideogramheight <- 10
# kp <- plotKaryotype(genome=genomeBuild, chromosomes=c("chr8"), plot.type=1, plot.params = pp)
# kpDataBackground(kp, r1=0.47)
# kpDataBackground(kp, r0=0.53)
# kpAddBaseNumbers(kp, cex=0.75)
# kpAxis(kp, r0=0.53, ymin=-4, ymax=4, numticks = 5, cex=1)
# kpAxis(kp, r1=0.47, ymin=-4, ymax=4, numticks = 5, cex=1)
# mtext(text="Off-Target Log Ratio", )
# kpPoints(kp, r0=0.53, chr=chr.ideo, x=data1.x, y=data1.y, col="blue", ymin=-4, ymax=4)
# kpPoints(kp, r1=0.47, chr=chr.ideo, x=data2.x, y=data2.y, col="red", ymin=-4, ymax=4)

}

save.image(outImage)
