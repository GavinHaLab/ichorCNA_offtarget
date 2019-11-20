# *Snakemake workflow for ichorCNA_offtarget*
Estimating tumor fraction by exploiting the off-target reads in targeted DNA sequencing.

# Description
The ichorCNA_offtarget.snakefile workflow will run the ichorCNA off-target pipeline, starting from bam files and generating files and plots documenting copy number alterations across the genome.

The plotCNAzoom_offtarget.snakefile workflow will use the results of ichorCNA_offtarget.snakefile and zoom in on particular regions of interest, producing plots that provide a closer look at copy number alterations in those regions.

# Requirements
* R-3.3 or above
  * ichorCNA
  * HMMcopy
  * optparse
* Python-3.4 or above
  * snakemake-3.12.0 or above
  * PyYAML-3.12 or above
* HMMcopy Suite (http://compbio.bccrc.ca/software/hmmcopy/)
  * In particular, readCounter is used.
* ichorCNA repository (https://github.com/broadinstitute/ichorCNA)
 
# Set-up
## config/samples.yaml
Please specify the samples to be analyzed in config/samples.yaml, following the format explained therein.
 
## config/config.yaml
There are a number of parameters to adjust in config/config.yaml.  Filepaths to where your ichorCNA repository (https://github.com/broadinstitute/ichorCNA) clone was downloaded to must be inserted in a few places, as well as the filepath to your readCounter binary.
A nice description of some of the other parameters that can be adjusted, such as bin size and segmentation settings, can be found in the README.md for the snakemake workflow of regular (not off-target) ichorCNA, at https://github.com/broadinstitute/ichorCNA/wiki/SnakeMake-pipeline-for-ichorCNA.

## formatting of .bed file containing targeted regions
Currently, this code handles .bed files of targeted regions with no headers and at least 5 columns.  Please ensure your .bed file follows this formatting, or change the following code in code/ichorCNA_offtarget.R to fit your purposes:

`baits <- fread(baitFile.tum, skip=0)`

`colnames(baits) <- c("chr", "start", "end", "name", "strands")`
 
# Running the pipeline
To invoke the snakemake workflows locally, use

`snakemake -s ichorCNA_offtarget.snakefile`

and

`snakemake -s plotCNAzoom_offtarget.snakefile`

To run the snakemake workflows on a slurm cluster, use

`snakemake -s ichorCNA_offtarget.snakefile --latency-wait 60 --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 50`

and

`snakemake -s plotCNAzoom_offtarget.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 50`

To perform a dry run, where the commands to be executed will be printed, but the snakefile will not actually be run, add `-np` at the end of these commands.  This helps to verify the workflow will run as expected.
