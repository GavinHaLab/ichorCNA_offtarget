# file:   ichorCNA.R
# author: Gavin Ha, Ph.D.
#         Fred Hutch
# contact: <gha@fredhutch.org>
#
# ichorCNA: https://github.com/broadinstitute/ichorCNA
# date:  July 24, 2019
# description: ichorCNA analysis of on- or off-target panel sequencing data. This script accepts as input already normalized log ratio data.

library(optparse)

option_list <- list(
  make_option(c("--logRFile"), type = "character", help = "Path to on- or off-target tsv file containing normalized log ratio values. Required."),
  make_option(c("--statsFile"), type = "character", default=NULL, help = "Path to tsv file containing on- and off-target statistics, including sex of individual. Default: [%default]"),
  make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
  make_option(c("--centromere"), type="character", default=NULL, help = "File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package. Default: [%default]"),
  make_option(c("--normalPanel"), type="character", default=NULL, help="Median corrected depth from panel of normals. Default: [%default]"),
  make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
  make_option(c("--normal"), type="character", default="0.5", help = "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [%default]"),
  make_option(c("--scStates"), type="character", default="NULL", help = "Subclonal states to consider. Default: [%default]"),
  make_option(c("--coverage"), type="numeric", default=NULL, help = "PICARD sequencing coverage. Default: [%default]"),
  make_option(c("--lambda"), type="character", default="NULL", help="Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data. Default: [%default]"),
  make_option(c("--lambdaScaleHyperParam"), type="numeric", default=3, help="Hyperparameter (scale) for Gamma prior on Student's-t precision. Default: [%default]"),
  #	make_option(c("--kappa"), type="character", default=50, help="Initial state distribution"),
  make_option(c("--ploidy"), type="character", default="2", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
  make_option(c("--maxCN"), type="numeric", default=7, help = "Total clonal CN states. Default: [%default]"),
  make_option(c("--estimateNormal"), type="logical", default=TRUE, help = "Estimate normal. Default: [%default]"),
  make_option(c("--estimateScPrevalence"), type="logical", default=TRUE, help = "Estimate subclonal prevalence. Default: [%default]"),
  make_option(c("--estimatePloidy"), type="logical", default=TRUE, help = "Estimate tumour ploidy. Default: [%default]"),
  make_option(c("--maxFracCNASubclone"), type="numeric", default=0.7, help="Exclude solutions with fraction of subclonal events greater than this value. Default: [%default]"),
  make_option(c("--maxFracGenomeSubclone"), type="numeric", default=0.5, help="Exclude solutions with subclonal genome fraction greater than this value. Default: [%default]"),
  make_option(c("--minSegmentBins"), type="numeric", default=50, help="Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
  make_option(c("--altFracThreshold"), type="numeric", default=0.05, help="Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [%default]"),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases. Default: [%default]"),
  make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params. Default: [%default]"),
  make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze. Default: [%default]"),
  make_option(c("--genomeBuild"), type="character", default="hg19", help="Geome build. Default: [%default]"),
  make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
  make_option(c("--normalizeMaleX"), type="logical", default=TRUE, help = "If male, then normalize chrX by median. Default: [%default]"),
  make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male. Default: [%default]"),
  make_option(c("--includeHOMD"), type="logical", default=FALSE, help="If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [%default]"),
  make_option(c("--txnE"), type="numeric", default=0.9999999, help = "Self-transition probability. Increase to decrease number of segments. Default: [%default]"),
  make_option(c("--txnStrength"), type="numeric", default=1e7, help = "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [%default]"),
  	make_option(c("--plotFileType"), type="character", default="pdf", help = "File format for output plots. Default: [%default]"),
	make_option(c("--plotYLim"), type="character", default="c(-4,6)", help = "ylim to use for chromosome plots. Default: [%default]"),
  make_option(c("--outDir"), type="character", default="./", help = "Output Directory. Default: [%default]"),
  make_option(c("--libdir"), type = "character", default=NULL, help = "Script library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files. Default: [%default]")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

library(HMMcopy)
library(GenomicRanges)
library(GenomeInfoDb)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

patientID <- opt$id
logRFile <- opt$logRFile
statsFile <- opt$statsFile
gcWig <- opt$gcWig
mapWig <- opt$mapWig
normal_panel <- opt$normalPanel
#exons.bed <- opt$exons.bed  # "0" if none specified
centromere <- opt$centromere
flankLength <- opt$rmCentromereFlankLength
normal <- eval(parse(text = opt$normal))
scStates <- eval(parse(text = opt$scStates))
lambda <- eval(parse(text = opt$lambda))
lambdaScaleHyperParam <- opt$lambdaScaleHyperParam
estimateNormal <- opt$estimateNormal
estimatePloidy <- opt$estimatePloidy
estimateScPrevalence <- opt$estimateScPrevalence
maxFracCNASubclone <- opt$maxFracCNASubclone
maxFracGenomeSubclone <- opt$maxFracGenomeSubclone
minSegmentBins <- opt$minSegmentBins
altFracThreshold <- opt$altFracThreshold
ploidy <- eval(parse(text = opt$ploidy))
coverage <- opt$coverage
maxCN <- opt$maxCN
txnE <- opt$txnE
txnStrength <- opt$txnStrength
normalizeMaleX <- as.logical(opt$normalizeMaleX)
includeHOMD <- as.logical(opt$includeHOMD)
fracReadsInChrYForMale <- opt$fracReadsInChrYForMale
chrXMedianForMale <- -0.1
outDir <- opt$outDir 
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
outPlotDir <- paste0(outDir, "/", patientID)
dir.create(outPlotDir, showWarnings = FALSE)
libdir <- opt$libdir
plotFileType <- opt$plotFileType
plotYLim <- eval(parse(text=opt$plotYLim))
gender <- NULL
outImage <- paste0(outDir,"/", patientID, ".RData")
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
chrs <- as.character(eval(parse(text = opt$chrs)))
chrTrain <- as.character(eval(parse(text=opt$chrTrain))); 
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize))); 
seqlevelsStyle(chrs) <- genomeStyle
seqlevelsStyle(chrNormalize) <- genomeStyle
seqlevelsStyle(chrTrain) <- genomeStyle

## load ichorCNA library or source R scripts
if (!is.null(libdir) && libdir != "None"){
	source(paste0(libdir,"/R/utils.R"))
	source(paste0(libdir,"/R/segmentation.R"))
	source(paste0(libdir,"/R/EM.R"))
	source(paste0(libdir,"/R/output.R"))
	source(paste0(libdir,"/R/plotting.R"))
} else {
    library(ichorCNA)
}

## load seqinfo 
seqinfo <- getSeqInfo(genomeBuild, genomeStyle)

save.image(outImage)

## LOAD IN ON AND OFF-TARGET GC-CORRECTED FILES ##
id <- patientID
if (substr(logRFile,nchar(logRFile)-2,nchar(logRFile)) == "txt") {
  logRFiles <- data.frame(cbind(patientID, logRFile))
} 

## LOAD OFF TARGET STATS FILE TO GET GENDER ##
stats <- read.delim(statsFile, header = TRUE, sep = "\t")
gender <- stats$Sex[1]

numSamples <- nrow(logRFiles)
logR.data <- GRangesList()
for (i in 1:numSamples) {
	id <- logRFiles[i, 1]
	message("Loading normalized log ratios: ", id)
	logR.in <- read.delim(logRFiles[i, 2], header = TRUE, sep = "\t")
	logR.in <- as(logR.in , "GRanges")
	logR.in$valid <- logR.in$valid.tumor
	logR.in$copy <- logR.in$copy.norm
	logR.in$reads <- logR.in$reads.tumor
	logR.in$cor.gc <- logR.in$cor.gc.tumor
	logR.in$cor.map <- logR.in$cor.map.tumor
	logR.data[[id]] <- logR.in
}

chrInd <- as.character(seqnames(logR.data[[1]])) %in% chrTrain
## get positions that are valid
valid <- logR.data[[1]]$valid
if (numSamples >= 2) {
  for (i in 2:length(logR.data)){ 
    valid <- valid & logR.data[[i]]$valid
  } 
}
save.image(outImage)

### RUN HMM ###
## store the results for different normal and ploidy solutions ##
ptmTotalSolutions <- proc.time() # start total timer
results <- list()
loglik <- as.data.frame(matrix(NA, nrow = length(normal) * length(ploidy), ncol = 7, 
                 dimnames = list(c(), c("init", "n_est", "phi_est", "BIC", 
                 												"Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))
counter <- 1
compNames <- rep(NA, nrow(loglik))
mainName <- rep(NA, length(normal) * length(ploidy))
#### restart for purity and ploidy values ####
for (n in normal){
  for (p in ploidy){
    if (n == 0.95 & p != 2) {
        next
    }
    logR <- as.data.frame(lapply(logR.data, function(x) { x$copy })) # NEED TO EXCLUDE CHR X #
    param <- getDefaultParameters(logR[valid & chrInd, , drop=F], maxCN = maxCN, includeHOMD = includeHOMD, 
                ct.sc=scStates, ploidy = floor(p), e=txnE, e.same = 50, strength=txnStrength)
    param$phi_0 <- rep(p, numSamples)
    param$n_0 <- rep(n, numSamples)
    
    ############################################
    ######## CUSTOM PARAMETER SETTINGS #########
    ############################################
    # 0.1x cfDNA #
    if (is.null(lambda)){
			logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
			param$lambda <- rep(logR.var, length(param$ct))
			param$lambda[param$ct %in% c(2)] <- logR.var 
			param$lambda[param$ct %in% c(1,3)] <- logR.var 
			param$lambda[param$ct >= 4] <- logR.var / 5
			param$lambda[param$ct == max(param$ct)] <- logR.var / 15
			param$lambda[param$ct.sc.status] <- logR.var / 5
    }else{
			param$lambda[param$ct %in% c(2)] <- lambda[2]
			param$lambda[param$ct %in% c(1)] <- lambda[1]
			param$lambda[param$ct %in% c(3)] <- lambda[3]
			param$lambda[param$ct >= 4] <- lambda[4]
			param$lambda[param$ct == max(param$ct)] <- lambda[2] / 15
			param$lambda[param$ct.sc.status] <- lambda[2] / 10
		}
		param$alphaLambda <- rep(lambdaScaleHyperParam, length(param$ct))  
    # 1x bulk tumors #
    #param$lambda[param$ct %in% c(2)] <- 2000
    #param$lambda[param$ct %in% c(1)] <- 1750
    #param$lambda[param$ct %in% c(3)] <- 1750
    #param$lambda[param$ct >= 4] <- 1500
    #param$lambda[param$ct == max(param$ct)] <- 1000 / 25
		#param$lambda[param$ct.sc.status] <- 1000 / 75
		#param$alphaLambda[param$ct.sc.status] <- 4
		#param$alphaLambda[param$ct %in% c(1,3)] <- 5
		#param$alphaLambda[param$ct %in% c(2)] <- 5
		#param$alphaLambda[param$ct == max(param$ct)] <- 4
				
		#############################################
		################ RUN HMM ####################
		#############################################
    hmmResults.cor <- HMMsegment(logR.data, valid, dataType = "copy", 
                                 param = param, chrTrain = chrTrain, maxiter = 50,
                                 estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                                 estimateSubclone = estimateScPrevalence, verbose = TRUE)
                                     
    for (s in 1:numSamples){
		iter <- hmmResults.cor$results$iter
		id <- names(hmmResults.cor$cna)[s]

		# correct integer copy number based on estimated purity and ploidy
   		correctedResults <- correctIntegerCN(cn = hmmResults.cor$cna[[s]],
   				segs = hmmResults.cor$results$segs[[s]], 
    			purity = 1 - hmmResults.cor$results$n[s, iter], ploidy = hmmResults.cor$results$phi[s, iter],
    			cellPrev = 1 - hmmResults.cor$results$sp[s, iter], 
    			maxCNtoCorrect.autosomes = maxCN, maxCNtoCorrect.X = maxCN, minPurityToCorrect = 0.03, 
    			gender = gender, chrs = chrs, correctHOMD = includeHOMD)
		hmmResults.cor$results$segs[[s]] <- correctedResults$segs
		hmmResults.cor$cna[[s]] <- correctedResults$cn
		## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
		## check if there is an altered segment that has at least a minimum # of bins
		segsS <- hmmResults.cor$results$segs[[s]]
		segsS <- segsS[segsS$chr %in% chrTrain, ]
		segAltInd <- which(segsS$event != "NEUT")
		maxBinLength = -Inf
		if (sum(segAltInd) > 0){
			maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
			maxSegRD <- GRanges(seqnames=segsS$chr[segAltInd[maxInd]], 
								ranges=IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
			hits <- findOverlaps(query=maxSegRD, subject=logR.data[[s]][valid, ])
			maxBinLength <- length(subjectHits(hits))
		}
		## check if there are proportion of total bins altered 
		# if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
		cnaS <- hmmResults.cor$cna[[s]]
		altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
		altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)
		if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)){
			hmmResults.cor$results$n[s, iter] <- 1.0
		}
      	## plot solution ##
		outPlotFile <- paste0(outPlotDir, "/", id, "_genomeWide_", "n", n, "-p", p)
		mainName[counter] <- paste0(id, ", n: ", n, ", p: ", p, ", log likelihood: ", signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4))
		plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType=plotFileType, 
					 plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, seqinfo=seqinfo, main=mainName[counter])
    }
    iter <- hmmResults.cor$results$iter
    results[[counter]] <- hmmResults.cor
    loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
    subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }))
    fracGenomeSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }))
    fracAltSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$copy.number != 2) }))
    fracAltSub <- lapply(fracAltSub, function(x){if (is.na(x)){0}else{x}})
    loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub, digits=2), collapse=",")
    loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSub), digits=2), collapse=",")
    loglik[counter, "init"] <- paste0("n", n, "-p", p)
    loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
    loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")

    counter <- counter + 1
  }
}
## get total time for all solutions ##
elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
message("Total ULP-WGS HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")

### SAVE R IMAGE ###
save.image(outImage)
#save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))

### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###
if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
	fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone & 
						 		   loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone)
	if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
		ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing=TRUE)]
	}else{ # otherwise just take largest likelihood
		ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
	}
}else{#sort by likelihood only
  ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
}

#new loop by order of solutions (ind)
outPlotFile <- paste0(outPlotDir, "/", id, "_genomeWide_all_sols")
for(i in 1:length(ind)) {
  hmmResults.cor <- results[[ind[i]]]
  turnDevOff <- FALSE
  turnDevOn <- FALSE
  if (i == 1){
  	turnDevOn <- TRUE
  }
  if (i == length(ind)){
  	turnDevOff <- TRUE
  }
  plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
                     plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, 
                     seqinfo = seqinfo,
                     turnDevOn = turnDevOn, turnDevOff = turnDevOff, main=mainName[ind[i]])
}

hmmResults.cor <- results[[ind[1]]]
hmmResults.cor$results$loglik <- as.data.frame(loglik)
hmmResults.cor$results$gender <- gender
hmmResults.cor$results$chrYCov <- NULL
hmmResults.cor$results$chrXMedian <- NULL
hmmResults.cor$results$coverage <- NULL

outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
                      results = hmmResults.cor$results, patientID = patientID, outDir=outDir)
outFile <- paste0(outDir, "/", patientID, ".params.txt")
outputParametersToFile(hmmResults.cor, file = outFile)

## plot solutions for all samples 
plotSolutions(hmmResults.cor, logR.data, chrs, outDir, numSamples=numSamples,
              plotFileType=plotFileType, plotYLim=plotYLim, 
              logR.column = "logR", call.column = "Corrected_Call",
              estimateScPrevalence=estimateScPrevalence, maxCN=maxCN)

