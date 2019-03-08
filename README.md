# YARP -- Yet Another RNA-seq analysis Pipeline

Author: [Fred P. Davis](http://fredpdavis.com)

## Contributors

- Alejandro Villarino contributed R code to perform functional analysis using clusterProfiler

## Overview

The goal of this package is to help you analyze your RNA-seq data in 
(1) well-documented, (2) easy to reproduce, and (3) standardized way.

## 1. What do I need to get started?

YARP is written for the [HPC computing environment at NIH](http://hpc.nih.gov).
You need to [request an account on the biowulf cluster](https://hpc.nih.gov/docs/accounts.html),
and sufficient disk space to store the raw sequence data and resulting analyses.
Disk space varies by how deeply your samples were sequenced, but ~5GB/sample is a
reasonable estimate.


## 2. General advice for computational work

### 2.1. Take responsibility.

Neither the experiment nor the analysis is complicated -- thousands of
people have done RNA-seq before -- but, you still have to invest time to
fully understand the experiment and how its measurements are interpreted.
If you can't convincingly explain what's going on in the experiment or
the analysis, keep reading or asking others until you can.
 
YARP is meant to help you get started -- not to offload your
responsibility, or to substitute for your thinking.

### 2.2. Write everything down.

Meticulous notes are just as critical for computational work as they are for
experimental work. I keep records in three ways:

1. Scripts. These are text files containing the actual commands that perform
each analysis step. Scripts force you to keep track of your analyses.
 
2. README files. If there is anything important to know or remember about
what happened in a specific directory, I write it in a README text file
that I keep in that specific directory.
 
3. Electronic lab notebook. I keep a chronological lab notebook made of
text file per day (see [mdlabbook](http://github.com/fredpdavis/mdlabbook)).
You should use whatever system you are comfortable with, but I highly
recommend you keep electronic notes. When you want to remember what commands
you tried, or what analysis you were working on Monday two weeks ago, you
are unlikely to remember accurately unless you can look back at notes.

### 2.3. Trust no one, especially yourself

This is of course an exaggeration -- ultimately you are trusting multiple
layers of software that operate the sequencer through to the programs
that generate the figures you will interpret. Be skeptical of all results.
Don't over-interpret or assume your results are correct: lots of things can
and do go wrong. Do you have positive and negative controls? If something
looks weird, it probably is.

### 2.4. Look inside.

Scripts are just text files -- view them in a text editor, or at the command
line using the `less FILENAME` command, to see the commands that are being
run.

### 2.5. Be tidy.

#### i. keep an organized directory structure.

I make a new folder for each project I work on, with the following structure:

- README -- file where I describe the project and broad goals
- data/ - directory storing raw data, each in its own folder -- eg, data/fastq/
- src/ - directory storing all scripts / code used in the project. I usually organize by language -- eg, src/R/
- run/ - directory where I run scripts, each kept in a dated directory -- eg, run/20181027.align_reads/
- results/ - store results from off-the-shelf programs -- eg, results/kallisto/
- analysis/ - where I do 'secondary' analysis, like making figures

YARP comes with, or builds up, the directory structure above.

#### ii. don't use weird characters

Don't use spaces or weird characters (asterisks, slashes, apostrophes,
quotes, etc) in your sample names, condition names, etc. Use underscores and
periods if you need separators.

#### iii. don't edit stuff unless you understand what you're doing.

#### iv. Turn off Mac autocorrect.

Apostrophes and quotes have specific meanings in scripts. By default, Macs will
often change these kinds of characters to curly quotes that will break the
command -- this is a problem when eg, you try to copy and paste a command into
a script that you are editing, or into the terminal to test a command.

## 3. RNA-seq analysis

The scripts in this package will perform the most routine (and often most informative) analyses:

1. estimating transcript abundance (= expression level)

2. identifying differentially expressed genes between pairs of conditions

3. visualizing results using genome browsers, scatterplots and heatmaps.

We will not cover the many other kinds of analyses you can perform on RNA-seq
measurements, including identifying alternative splicing events, estimating
nascent transcription, or evaluating more complex experimental designs.

YARP uses the following external data sources. __You don't need to
download these__, as this package either comes with or will retrieve the
appropriate files.

| data set                                                                                      | purpose                                       |
|-----------------------------------------------------------------------------------------------|-----------------------------------------------|
| [ENSEMBL GRCm38.94 release94](https://useast.ensembl.org/Mus_musculus/Info/Index)             | genome sequence, gene annotations             |
| [Life Tech: ERCC92](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip)        | synthetic spike-in sequences, concentratsion  |


YARP uses the following underlying tools to analyze RNA-seq data.
__You don't need to download these__, as they are already installed on the NIH
HPC through their fantastic system of ['modules'](https://hpc.nih.gov/apps/modules.html)

| software                                                      | purpose                                                   |
|---------------------------------------------------------------|-----------------------------------------------------------|
| [kallisto 0.44.0](https://pachterlab.github.io/kallisto/)     | pseudoalign reads to transcriptome, estimates abundance   |
| [STAR 2.5.4a](https://github.com/alexdobin/STAR)              | align reads to genome                                     |
| [deeptools 3.1.2](https://deeptools.readthedocs.io/)          | make bigwig files for visualizing in IGV                  |
| [sleuth 0.29](https://pachterlab.github.io/sleuth/)           | R package for differential expression analysis            |
| [R 3.5.0](https://www.r-project.org)                          | creating figures, and running sleuth                      |

These programs offer lots of options -- this package uses particular
incantations of each, but as you get more comfortable with the analysis, look
inside the scripts to see what specific options are used for each program, and
then read the program manuals to understand what each does and what other
options you may want to try.


### 3.1. Setup your project directory

- Login to biowulf

    - on Mac: open the Terminal program, and type `ssh biowulf.nih.gov`, press enter, and enter your password when prompted

    - on Windows: download and run [PuTTY](https://www.putty.org) to ssh into biowulf.nih.gov, specifying your username and password

- if you don't already have one, make a projects directory on your /data share.

    - NOTE: __REPLACE davisfp with your username__

```
mkdir /data/davisfp/projects
```

- change directories to your projects directory

```
cd /data/davisfp/projects
```

- download this package

```
git clone https://github.com/fredpdavis/oshealab_genomic_scripts.git
```

- rename to your project name, cytokineX in this example. __NOTE: REPLACE
cytokineX with your project name__

```
mv oshealab_genomic_scripts cytokineX
```

- change directories to this project.

```
cd cytokineX
```

- List the files to see what you're getting in this package.

```
tree .
```

### 3.2. Prepare your input files

Now that you have a set of scripts and basic directory structure, you can start
using it to analyze your own data.

#### 3.2.1 Get FASTQ sequence files

For testing, I provide four small fastq files -- you can delete these if you want.

```
shell> ls data/fastq/test
    
in1.fq.gz
in2.fq.gz
in3.fq.gz
in4.fq.gz
```

Sequence files often have a .fastq or .fq suffix, and are typically
compressed by either gzip or bunzip, resulting in a further '.gz' or
'.bz2' suffix, respectively.

You can request your own sequence files from the NIAMS core if your samples were
sequenced locally, or download them from
[NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) if you want to re-process
published data.

Store the sequence files in the `data/fastq/` directory in subdirectories
named by sequencing run identifiers (if sequenced locally), or by author if
the FASTQ comes from a published paper.


#### 3.2.2 Create the sample sheet `metadata/rnaseq_samples.txt`

- This file describes all of your RNA-seq samples.

- This is a tab-delimited text file. Make sure you use tab's, __NOT SPACES__

- Each line represents an RNA-seq sample.

- This package includes an example file listing the test samples

```
shell> cat metadata/test_rnaseq_samples.txt

sampleID	sampleName	runID	fastqName	cytokineStim
r001	unstim_rep1	test	in1.fq.gz	nc
r002	unstim_rep3	test	in2.fq.gz	nc
r003	IFNg_72h_rep1	test	in3.fq.gz	IFNg_72h
r004	IFNg_72_rep2	test	in4.fq.gz	IFNg_72h
```

- Make a copy of this file and edit it to describe your samples. You can do the latter using nano, an easy to use text editor.

```
cp metadata/test_rnaseq_samples.txt metadata/rnaseq_samples.txt
nano metadata/rnaseq_samples.txt
```

- requires 4 fields:
    1. sampleID
    2. sampleName
    3. runID
    4. fastqName

- add additional columns to describe sample properties like celltype,
  tissue, etc. this is useful for specifying the condition pairs you
  want to compare in the next file.

- expects to find fastq files in `data/fastq/<runID>/<fastqName>`.
For example, the read data for sampleID r001 should be in `data/fastq/test/in1.fq.gz`.


#### 3.2.3 Create the comparisons file `metadata/rnaseq_comparisons.txt`

- This file lists groups of samples that you'd like to compare

- This is a 4-column tab-delimited text file. Make sure you use tab's, __NOT SPACES__

- Each line represents a pairwise comparison

- Samples can be specified by sampleName or by other features specified in
the sample sheet

```
shell> cat metadata/test_rnaseq_comparisons.txt

group1name      group1criteria  group2name      group2criteria
no.cytokine     cytokineStim=nc gamma   cytokineStim=IFNg_72h
```

- Make a copy of this file and edit it to describe the comparisons you'd like to make:

```
cp metadata/test_rnaseq_comparisons.txt metadata/rnaseq_comparisons.txt
nano metadata/rnaseq_samples.txt
```

- If you want to specify your samples using more than one feature, use
commas to express logical AND, and semicolon to express logical OR.
AND's will be interpreted first. For example, `cytokineStim=nc,cellType=CD4;celltype=CD8`
will pick either CD4 samples without cytokine stimulation OR any CD8 samples.

### 4. Retrieve and prepare genome and gene information

This step will retrieve the files necessary to process mouse data. Take a look
inside the script to get a sense of what data is being retrieved. If you ever
want to change genome version, gene annotations, species, etc. you will edit
this file.

- By default, YARP uses ENSEMBL release 94. If you want to eventually compare
your bulk RNA-seq to 10x single cell RNA-seq that you process with cellranger,
you should instead use the sambe ENSEMBL version. Cellranger version 1, 2, and 3
use ENSEMBL release 84, 84, and 93, respectively.

- If necessary, you can specify the ENSEMBL release by editing the following
three files at the specified lines. Simply change 94 to your desired ENSEMBL
release.

    - `src/r/basicRnaSeqAnalysis.R` line 498: `specs$ensembl_release<- "94"`
    - `src/slurm/prepare_indices.slurm.csh` line 21 `set ENSEMBL_RELEASE="94"`
    - `src/slurm/process_rnaseq_samples.slurm.csh` line 43 `set ENSEMBL_RELEASE="94"`

- Make a new directory and copy the prepared_indices.slurm.csh script there.

```
mkdir -p run/20181102.prepare_inputs
cp src/slurm/prepare_indices.slurm.csh run/20181102.prepare_inputs
```

- Switch to that directory

```
cd run/20181102.prepare_inputs
```

- Edit the script (eg with nano: `nano prepare_indices.slurm.csh`) to specify the BASEDIR. That is, change cytokineX to the name of your project:

```
set BASEDIR="/data/davisfp/projects/cytokineX"
```

- Make a directory to hold the output of the script

```
mkdir slurm_out
```

- Submit the script to the cluster to download and process necessary files:

```
sbatch prepare_indices.slurm.csh
```

- To check if the job is complete, use the squeue command

```
squeue -u davisfp
```

You only need to run this step once (per project). You don't need to run this
step again if you just want to process additional samples.

### 5. Process your samples

- Make a new run directory and copy the next script there.

```
cd ..
mkdir 20181102.process_samples
cp ../src/slurm/process_rnaseq_samples.slurm.csh 20181102.process_samples
```

- Edit the script (eg with nano: `nano process_rnaseq_samples.slurm.csh`) to specify the BASEDIR. That is, change cytokineX to the name of your project:

```
set BASEDIR="/data/davisfp/projects/cytokineX"
```

- Edit the script to specify which tasks to process. By default, all samples
listed in the metadata file will be processed. If you just want to process
a subset of those samples, you can specify their defining features in the
script. For example, to process only samples not stimulated by cytokine:

```
set SAMPLE_OPTION="-cytokineStim no"
```

- Edit the script to tell the cluster how many jobs you will run. For example,
  we'd like to process the 4 samples in our test set:

```
#SBATCH --array=1-4
```

- It makes sense to run this script on a single job to begin with, just to
  make sure that everything runs properly:

```
#SBATCH --array=1-1
```

- Make a directory to hold the script output:

```
mkdir slurm_out
```

- Submit the jobs to the cluster:

```
sbatch process_rnaseq_samples.slurm.csh
```

- Check job status with `squeue`

- Browse the results directory to see the files that were generated:

### 6. Create figures and tables

Now that the samples have been individually processed, we will next make some
standard figures to see how the samples relate to one another, identify
differentially expresed genes and create figures.

This step is implemented in an [R](https://www.r-project.org/) script. Unlike
the slurm shell scripts in the previous steps, which we copied into the run
directory to edit, we won't need to edit the R script and so will just run it
where it sits in `src/R/basicRnaSeqAnalysis.R`.

We will run R on the cluster -- instead of submitting a "batch" job with sbatch,
we will use the sinteractive command to request an interactive session, meaning
a command line prompt on a computer in the cluster.

- make an analysis directory and switch to it.

```
cd ../../
mkdir -p analysis/20181102.makeFigures
cd analysis/20181102.makeFigures
```

- request an interactive session

```
sinteractive --x11 --cpus-per-task=2 --mem=64g --time=24:00:00
```

- load the R 'module'

```
module load R
```

- start R 

```
R
```

- Load the R script

```
> source("../../src/R/basicRnaSeqAnalysis.R")
```

- Load the raw data, by running the main() routine, and providing the location of your base directory:

```
> dat <- main(baseDir = "~/data/projects/cytokineX", returnData=TRUE)
```

- If you're curious what you just loaded, you can use the str() command to see what is in dat:


```
> str(dat)
```

- Load additional data ( differential expression calls)

```
> dat <- main(dat, returnData=TRUE)
```

- Make a basic set of figures:

```
> makeFigures(dat)
```

- This command will make these basic figures:

    1. heatmap of variably expressed genes across all samples. hierarchicaly clustered.
    2. heatmap of correlation between transcript abundances across all pairs of samples.
    3. scatterplots of estimated transcript abundance between pairs of samples specified in rnaseq_comparisons.txt, highlighting differentially expressed genes
    4. heatmap of differentially expressed genes

- You can also customize the figures to make them more useful. For example, if you'd like to annotate the correlation heatmap with sample properties (specified in the rnaseq_samples.txt file):

```
makeCorrHeatmap(dat, sampleAnnotate=c("cytokineStim"))

```

This will make a heatmap where the rows and columns are annotated with different colors if they were stimulated with cytokine or not.

- By default, the plotHmap.deg() routine will show all differentially expressed genes, but you can also specify your own genes like this:

```
plotHmap.deg(dat, geneList = c("Jak1", "Jak2", "Jak3", "Jak4"))
```

- You can also annotate the plotHmap.deg() heatmap with sample properties, just like for the correlation heatmap

```
plotHmap.deg(dat, sampleAnnotate=c("cytokineStim"))
```

- By default, the plotHmap.deg() routine will show individual samples, but you can also show averages across groups of samples (eg, if you want to show replicate averages).

```
plotHmap.deg(dat, repAvg = list(unstim = c("unstim_rep1", "unstim_rep2"), gamma=c("IFNg_72h_rep1", "IFNg_72h_rep2")))
```

- To add sample properties to the variable gene heatmap (made by makeVarGeneHeatmap()), you can specify what properties to show:

```
makeVarGeneHeatmap(dat,sampleAnnotate=c("cytokineStim"))
```

- To view the figures, secury copy (scp) the whole directory to your local machine's Desktop. Open a new terminal window on your local machine and type:

```
scp -r biowulf:data/projects/cytokineX/analysis/basicFigures ~/Desktop
```

## Design decisions

### Options

- genome version. default mm10
- annotation source/version. default ENSEMBL release 94
- ERCC spike-ins. by default added to index. if not used in a sample, negligible contribution anyways and can always delete them.
- gene types. default only cDNA (protein-coding genes)


