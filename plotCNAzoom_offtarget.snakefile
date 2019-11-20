#plotCNAzoom_offtarget.snakefile
#Ha Lab
#Fred Hutchinson Cancer Research Center

configfile: "config/configPlotZoom_offtarget.yaml"
configfile: "config/samples.yaml"

import glob
def getTITANpath(base, id, ext):
  return glob.glob(''.join([base, "results/titan/optimalClusterSolution/", id, "_cluster*", ext]))


rule all:
  input:
  	expand("results/plotCNAzoom/{plotID}/{tumor}_CNA-SV_chr{chr}-{start}-{end}.{format}", tumor=config["pairings"], plotID=config["plot_id"], chr=config["plot_chr"], start=config["plot_startPos"], end=config["plot_endPos"], format=config["plot_format"])


rule plotCNAzoom:
	input:
		offCNA="results/ichorCNA/{tumor}/offTarget/{tumor}.cna.seg",
		offParams="results/ichorCNA/{tumor}/offTarget/{tumor}.params.txt",
		onCNA="results/ichorCNA/{tumor}/onTarget/{tumor}.cna.seg"
	output:
		"results/plotCNAzoom/{plotID}/{tumor}_CNA-SV_chr{chr}-{start}-{end}.{format}"
	params:
		plotCNscript=config["plotCN_script"],
		plotfuncs=config["plot_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		zoom=config["plot_zoom"],
		chr=config["plot_chr"],
		start=config["plot_startPos"],
		end=config["plot_endPos"],
		ylim=config["plot_ylim"],
		geneFile=config["plot_geneFile"],
		size=config["plot_size"],
		format=config["plot_format"]
	log:
		"logs/plotCNAzoom/{plotID}/{tumor}_chr{chr}-{start}-{end}.{format}.log"
	shell:
		"Rscript {params.plotCNscript} --id {wildcards.tumor} --plot_funcs {params.plotfuncs} --offTargetCNFile {input.offCNA} --onTargetCNFile {input.onCNA} --OffTargetParamFile {input.offParams} --chrs \"{params.chr}\" --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --start {params.start} --end {params.end} --zoom {params.zoom} --plotYlim \"{params.ylim}\" --geneFile {params.geneFile} --plotSize \"{params.size}\" --outPlotFile {output} > {log} 2> {log}"
