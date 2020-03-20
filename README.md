# workflows
Nextflow pipelines for performing bioinformatics workflows on lewis

## Prerequisites
### Installing nextflow
Nextflow is easy to install:
```bash
wget -qO- https://get.nextflow.io | bash
```
This command creates the `nextflow` binary in the current directory. You can
put it somewhere in your `$PATH` if you want.

### Configuring nextflow
Nextflow looks for a file called `nextflow.config` in the directory you are
running it from to figure out how you want to run it. Here is my configuration
file for running pipelines on lewis using SLURM:
```groovy
process {
    executor = 'slurm'
    queue = 'BioCompute'
    time = '2d'
    memory = '48 GB'
    clusterOptions = '--account=warrenlab'
}

singularity {
    enabled = true
    runOptions = "--bind /storage"
}
```
The `process` block tells nextflow to submit jobs to the cluster using SLURM
rather than trying to run them on the login node (important!), and to ask SLURM
for the listed queue, time, memory, and accounts. You can edit these for the
workflow you are running, and override any of them for a specific job if, for
example, there is one step that requires more memory than the others.

The `singularity` block tells nextflow to use
[Singularity](https://sylabs.io/singularity/) to run any processes that ask for
a container. Singularity is installed on lewis, so you don't need to do anything
special to run containerized processes besides add these lines to your config
file.

#### Some lewis-specific notes
Nextflow needs to be run from a directory that allows file locking. On lewis,
`htc` allows file locking but `hpc` doesn't. However, `hpc` is much faster, and
so it is where you should be doing all heavy I/O work. I get around this by:
1. having a project directory on `htc` where I keep the pipeline and config file
and run nextflow from
2. keeping all input data on `hpc`, and
3. setting the environment variable `NXF_WORK` to a directory on `hpc` where
nextflow will write all files

Then, at the end of the pipeline execution, nextflow will make links in your
`htc` directory to the output files on `hpc`.

Another thing: all software used here is either already installed on lewis, in
which case nextflow takes care of loading the necessary modules, or
containerized, in which case nextflow takes care of downloading the container
image, building the container, and running the commands inside the container, or
the software is on bioconda, in which case nextflow will make a new conda
environment with the necessary packages installed and run commands there. You
do need to have [anaconda](https://www.anaconda.com/) installed for that.

## Generic instructions on how to run a workflow
* Copy the workflow file (ending with `.nf`) to the directory where you'd like
to run nextflow from.
* Write a config file as described above and save it in the same directory as
`nextflow.config`.
* Edit the workflow file so that the parameters point to your input files. I
keep all parameters at the very beginning of each workflow so that you don't
need to dig through the file to set them.
* (optional) Edit the workflow file a little bit more to fine-tune parameters
like number of CPUs, amount of memory, versions of software installed on lewis,
etc.
* Run the pipeline with `nextflow run workflow.nf`. If the workflow uses
singularity, you'll have to do this from an interactive session as singularity
is disabled on the login node.

## Workflow descriptions and parameters
### arrow
Uses the [arrow](https://github.com/PacificBiosciences/gcpp) algorithm to
polish a PacBio assembly with long reads.

__Input:__
* an unpolished long-read genome assembly in fasta format
* a directory full of PacBio long reads in bam format, with the associated pbi
index files

__Output:__ a polished genome assembly in fasta format, in the directory
`arrow_out`

__Steps:__
* creates a PacBio-friendly reference dataset
* aligns all of the long reads to this reference using `pbmm2`
* splits the genome and alignments into chunks to paralellize
* runs `gcpp` in arrow mode to polish each chunk
* puts all the chunks back together

__Parameters:__
* `reference`: path to the genome assembly fasta file
* `subreads_folder`: path to a folder containing bam files and pbi files
* `numChunks`: number of chunks to divide the genome up into

__Notes:__
* it is recommended to polish a genome assembly twice. So, once you are done
running this pipeline the first time, take the output polished fasta file and
use it as input to a second iteration of this pipeline.

### manta
Uses the [manta](https://github.com/Illumina/manta) variant caller to call
structural variants using short reads.

__Input:__
* a genome assembly in fasta format
* paired-end short reads, either as a wildcard-expandable path to fastq files,
or a list of SRA accessions

__Output:__
* a vcf file with structural variants called for all samples, in `manta_out`
* a bwa index of your reference genome, in `bwa_index`
* an indexed bam file for each sample, in `alignments`

__Steps:__
* indexes the reference genome for bwa
* aligns all reads to the reference genome using bwa
* runs manta to call variants

__Parameters:__
* `reference`: path to the reference genome in fasta format
* reads: either use the line starting with `.fromFilePairs`, edit the string to
point to the path of all of your short read files, and leave the line starting
with `.fromSRA` commented out, __OR__ uncomment the line starting with
`.fromSRA`, edit it to point to a file containing a list of SRA accessions, one
per line, and comment out the line starting with `.fromFilePairs`.

### smoove
Uses the [smoove](https://github.com/brentp/smoove) wrapper for the
[lumpy](https://github.com/arq5x/lumpy-sv) variant caller to call structural
variants using short reads.

__Input:__
* a genome assembly in fasta format
* paired-end short reads, either as a wildcard-expandable path to fastq files,
or a list of SRA accessions

__Output:__
* a vcf file with structural variants called for all samples, in `output`
* a bwa index of your reference genome, in `bwa_index`
* an indexed bam file for each sample, in `alignments`

__Steps:__
* indexes the reference genome for bwa
* aligns all reads to the reference genome using bwa
* runs `smoove call` to call variants on each sample separately
* runs `smoove merge` to make a list of all the variants called in all samples
* runs `smoove genotype` to re-call variants for each sample, but jointly with
all the other samples this time
* runs `smoove paste` to put all the samples into a single vcf

__Parameters:__
* `reference`: path to the reference genome in fasta format
* reads: either use the line starting with `.fromFilePairs`, edit the string to
point to the path of all of your short read files, and leave the line starting
with `.fromSRA` commented out, __OR__ uncomment the line starting with
`.fromSRA`, edit it to point to a file containing a list of SRA accessions, one
per line, and comment out the line starting with `.fromFilePairs`.
* `scratch`: a directory for temporary files, preferably local to the compute
node rather than shared for speed purposes. On `BioCompute`, each node has a
scratch directory for each user at `/local/scratch/username`.


### trio-binning
TBD

### trio-scaffolding
TBD
