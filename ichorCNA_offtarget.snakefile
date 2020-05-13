configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
  input:
  	expand("results/readDepth/{samples}.bin{binSize}.wig", samples=config["samples"], binSize=str(config["binSize"])),
  	expand("results/normalizeOfftarget/{tumor}/{tumor}_offTarget_cor.txt", tumor=config["pairings"]),
	expand("results/normalizeOfftarget/{tumor}/{tumor}_onTarget_cor.txt", tumor=config["pairings"]),
	expand("results/normalizeOfftarget/{tumor}/{tumor}_allTarget_cor.txt", tumor=config["pairings"]),
	expand("results/normalizeOfftarget/{tumor}/{tumor}_readStats.txt", tumor=config["pairings"]),
	expand("results/ichorCNA/{tumor}/offTarget/{tumor}.cna.seg", tumor=config["pairings"]),
	expand("results/ichorCNA/{tumor}/onTarget/{tumor}.cna.seg", tumor=config["pairings"]),
	expand("results/ichorCNA/{tumor}/allTarget/{tumor}.cna.seg", tumor=config["pairings"])

rule read_counter:
	input:
		lambda wildcards: config["samples"][wildcards.samples]
	output:
		"results/readDepth/{samples}.bin{binSize}.wig"
	params:
		readCounter=config["readCounterScript"],
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"]
	resources:
		mem=4
	log:
		"logs/readDepth/{samples}.bin{binSize}.log"
	shell:
		"{params.readCounter} {input} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"


rule normalizeOffTarget:
	input:
		tum="results/readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
		norm=lambda wildcards: "results/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
	output:
		offNorm="results/normalizeOfftarget/{tumor}/{tumor}_offTarget_cor.txt",
		onNorm="results/normalizeOfftarget/{tumor}/{tumor}_onTarget_cor.txt",
		allNorm="results/normalizeOfftarget/{tumor}/{tumor}_allTarget_cor.txt",
		stats="results/normalizeOfftarget/{tumor}/{tumor}_readStats.txt"
	params:
		outDir="results/normalizeOfftarget/{tumor}/",
		rscript=config["offTarget_script"],
		offTargetFuncs=config["offTarget_utils_script"],
		id="{tumor}",
		chrs=config["ichorCNA_chrs"],
		baitBedTum=config["ichorCNA_targets"],
		gcWig=config["ichorCNA_gcWig"],
		mapWig=config["ichorCNA_mapWig"],
		mapScoreThres=config["ichorCNA_mapScoreThres"],
		centromere=config["ichorCNA_centromere"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"],
		libdir=config["ichorCNA_libdir"]
	log:
		"logs/normalizeOfftarget/{tumor}.log"
	shell:
		"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --offTargetFuncs {params.offTargetFuncs} --TUMWIG {input.tum} --NORMWIG {input.norm} --baitBedTum {params.baitBedTum} --gcWig {params.gcWig} --mapWig {params.mapWig} --mapScoreThres {params.mapScoreThres} --centromere {params.centromere} --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --fracReadsInChrYForMale {params.fracReadsChrYMale} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"

rule ichorCNA_OffTarget:
	input:
		logR="results/normalizeOfftarget/{tumor}/{tumor}_offTarget_cor.txt",
		stats="results/normalizeOfftarget/{tumor}/{tumor}_readStats.txt"
	output:
		#corrDepth="results/ichorCNA/{tumor}/offTarget/{tumor}.correctedDepth.txt",
		#param="results/ichorCNA/{tumor}/offTarget/{tumor}.params.txt",
		cna="results/ichorCNA/{tumor}/offTarget/{tumor}.cna.seg",
		#segTxt="results/ichorCNA/{tumor}/offTarget/{tumor}.seg.txt",
		#seg="results/ichorCNA/{tumor}/offTarget/{tumor}.seg",
		#rdata="results/ichorCNA/{tumor}/offTarget/{tumor}.RData"
	params:
		outDir="results/ichorCNA/{tumor}/offTarget/",
		rscript=config["ichorCNA_rscript"],
		id="{tumor}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
				sex=config["ichorCNA_sex"],
		gcwig=config["ichorCNA_gcWig"],
		mapwig=config["ichorCNA_mapWig"],
		repTimeWig=config["ichorCNA_repTimeWig"],
		normalpanel=config["ichorCNA_normalPanel"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		likModel=config["ichorCNA_likModel"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		centromere=config["ichorCNA_centromere"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		minMapScore=config["ichorCNA_minMapScore"],
		maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
		normal2IgnoreSC=config["ichorCNA_normal2IgnoreSC"],
		scPenalty=config["ichorCNA_scPenalty"],
		exons=config["ichorCNA_exons"],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"],
		libdir=config["ichorCNA_libdir"]
	resources:
		mem=4
	log:
		"logs/ichorCNA/offTarget/{tumor}.log"
	shell:
		"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --logRFile {input.logR} --statsFile {input.stats}  --gcWig {params.gcwig} --mapWig {params.mapwig} --repTimeWig {params.repTimeWig} --sex {params.sex} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --likModel {params.likModel} --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --normal2IgnoreSC {params.normal2IgnoreSC} --scPenalty {params.scPenalty} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"


rule ichorCNA_OnTarget:
	input:
		logR="results/normalizeOfftarget/{tumor}/{tumor}_onTarget_cor.txt",
		stats="results/normalizeOfftarget/{tumor}/{tumor}_readStats.txt"
	output:
		#corrDepth="results/ichorCNA/{tumor}/onTarget/{tumor}.correctedDepth.txt",
		#param="results/ichorCNA/{tumor}/onTarget/{tumor}.params.txt",
		cna="results/ichorCNA/{tumor}/onTarget/{tumor}.cna.seg",
		#segTxt="results/ichorCNA/{tumor}/onTarget/{tumor}.seg.txt",
		#seg="results/ichorCNA/{tumor}/onTarget/{tumor}.seg",
		#rdata="results/ichorCNA/{tumor}/onTarget/{tumor}.RData"
	params:
		outDir="results/ichorCNA/{tumor}/onTarget/",
		rscript=config["ichorCNA_rscript"],
		id="{tumor}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
				sex=config["ichorCNA_sex"],
		gcwig=config["ichorCNA_gcWig"],
		mapwig=config["ichorCNA_mapWig"],
		repTimeWig=config["ichorCNA_repTimeWig"],
		normalpanel=config["ichorCNA_normalPanel"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		likModel=config["ichorCNA_likModel"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		centromere=config["ichorCNA_centromere"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		minMapScore=config["ichorCNA_minMapScore"],
		maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
		normal2IgnoreSC=config["ichorCNA_normal2IgnoreSC"],
		scPenalty=config["ichorCNA_scPenalty"],
		exons=config["ichorCNA_exons"],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"],
		libdir=config["ichorCNA_libdir"]
	resources:
		mem=4
	log:
		"logs/ichorCNA/onTarget/{tumor}.log"
	shell:
		"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --logRFile {input.logR} --statsFile {input.stats}  --gcWig {params.gcwig} --mapWig {params.mapwig} --repTimeWig {params.repTimeWig} --sex {params.sex} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --likModel {params.likModel} --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --normal2IgnoreSC {params.normal2IgnoreSC} --scPenalty {params.scPenalty} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"


rule ichorCNA_AllTarget:
	input:
		logR="results/normalizeOfftarget/{tumor}/{tumor}_allTarget_cor.txt",
		stats="results/normalizeOfftarget/{tumor}/{tumor}_readStats.txt"
	output:
		#corrDepth="results/ichorCNA/{tumor}/onTarget/{tumor}.correctedDepth.txt",
		#param="results/ichorCNA/{tumor}/onTarget/{tumor}.params.txt",
		cna="results/ichorCNA/{tumor}/allTarget/{tumor}.cna.seg",
		#segTxt="results/ichorCNA/{tumor}/onTarget/{tumor}.seg.txt",
		#seg="results/ichorCNA/{tumor}/onTarget/{tumor}.seg",
		#rdata="results/ichorCNA/{tumor}/onTarget/{tumor}.RData"
	params:
		outDir="results/ichorCNA/{tumor}/allTarget/",
		rscript=config["ichorCNA_rscript"],
		id="{tumor}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
				sex=config["ichorCNA_sex"],
		gcwig=config["ichorCNA_gcWig"],
		mapwig=config["ichorCNA_mapWig"],
		repTimeWig=config["ichorCNA_repTimeWig"],
		normalpanel=config["ichorCNA_normalPanel"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		likModel=config["ichorCNA_likModel"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		centromere=config["ichorCNA_centromere"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		minMapScore=config["ichorCNA_minMapScore"],
		maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
		normal2IgnoreSC=config["ichorCNA_normal2IgnoreSC"],
		scPenalty=config["ichorCNA_scPenalty"],
		exons=config["ichorCNA_exons"],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"],
		libdir=config["ichorCNA_libdir"]
	resources:
		mem=4
	log:
		"logs/ichorCNA/allTarget/{tumor}.log"
	shell:
		"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --logRFile {input.logR} --statsFile {input.stats}  --gcWig {params.gcwig} --mapWig {params.mapwig} --repTimeWig {params.repTimeWig} --sex {params.sex} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --likModel {params.likModel} --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --normal2IgnoreSC {params.normal2IgnoreSC} --scPenalty {params.scPenalty} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"

