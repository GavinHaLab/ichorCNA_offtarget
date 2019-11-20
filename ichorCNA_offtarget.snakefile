#ichorCNA_offtarget.snakefile
#Ha Lab
#Fred Hutchinson Cancer Research Center

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
  input:
  	expand("results/readDepth/{samples}.bin{binSize}.wig", samples=config["samples"], binSize=str(config["binSize"])),
  	expand("results/normalizeOfftarget/{tumor}/{tumor}_offTarget_cor.txt", tumor=config["pairings"]),
	expand("results/normalizeOfftarget/{tumor}/{tumor}_onTarget_cor.txt", tumor=config["pairings"]),
	expand("results/normalizeOfftarget/{tumor}/{tumor}_readStats.txt", tumor=config["pairings"]),
	expand("results/ichorCNA/{tumor}/offTarget/{tumor}.cna.seg", tumor=config["pairings"]),
	expand("results/ichorCNA/{tumor}/onTarget/{tumor}.cna.seg", tumor=config["pairings"])

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
		stats="results/normalizeOfftarget/{tumor}/{tumor}_readStats.txt"
	params:
		outDir="results/normalizeOfftarget/{tumor}/",
		rscript=config["offTarget_script"],
		id="{tumor}",
		chrs=config["ichorCNA_chrs"],
		baitBedTum=config["ichorCNA_targets"],
		gcWig=config["ichorCNA_gcWig"],
		mapWig=config["ichorCNA_mapWig"],
		centromere=config["ichorCNA_centromere"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"],
		libdir=config["ichorCNA_libdir"]
	log:
		"logs/normalizeOfftarget/{tumor}.log"
	shell:
		"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --TUMWIG {input.tum} --NORMWIG {input.norm} --baitBedTum {params.baitBedTum} --gcWig {params.gcWig} --mapWig {params.mapWig} --centromere {params.centromere} --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"

rule ichorCNA_OffTarget:
	input:
		logR="results/normalizeOfftarget/{tumor}/{tumor}_offTarget_cor.txt",
		stats="results/normalizeOfftarget/{tumor}/{tumor}_readStats.txt"
	output:
		cna="results/ichorCNA/{tumor}/offTarget/{tumor}.cna.seg",
	params:
		outDir="results/ichorCNA/{tumor}/offTarget/",
		rscript=config["ichorCNA_rscript"],
		id="{tumor}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		normalpanel=config["ichorCNA_normalPanel"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		centromere=config["ichorCNA_centromere"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		maxFracGenomeSubclone=["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=["ichorCNA_maxFracCNASubclone"],
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
		"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --logRFile {input.logR} --statsFile {input.stats} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --txnE {params.txnE} --txnStrength {params.txnStrength} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"


rule ichorCNA_OnTarget:
	input:
		logR="results/normalizeOfftarget/{tumor}/{tumor}_onTarget_cor.txt",
		stats="results/normalizeOfftarget/{tumor}/{tumor}_readStats.txt"
	output:
		cna="results/ichorCNA/{tumor}/onTarget/{tumor}.cna.seg",
	params:
		outDir="results/ichorCNA/{tumor}/onTarget/",
		rscript=config["ichorCNA_rscript"],
		id="{tumor}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		normalpanel=config["ichorCNA_normalPanel"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		centromere=config["ichorCNA_centromere"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
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
		"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --logRFile {input.logR} --statsFile {input.stats} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --txnE {params.txnE} --txnStrength {params.txnStrength} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"
