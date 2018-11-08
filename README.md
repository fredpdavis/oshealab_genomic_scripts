# oshealab_genomic_scripts -- scripted and reproducible genomic data analysis

Any questions, contact fred.davis@nih.gov

## Overview

The goal of this package is to help you analyze your genomic data in 
(1) well-documented, (2) easy to reproduce, and (3) standard way 

## Getting started

This package is written for the HPC computing environment at NIH.
You need accounts on helix and the biowulf cluster, and sufficient disk space to
deal with sequence data and the resulting analyses. Disk space varies by how
deeply your samples were sequenced, but ~5GB/sample is a reasonable estimate.

Request your account [here](https://hpc.nih.gov/docs/accounts.html)

## General advice for computational work

### 1. Take responsibility.

Neither the experiment nor the analysis is complicated -- thousands of
people have done RNA-seq before -- but, you still have to invest time to
fully understand the experiment and how its measurements are interpreted.
If you can't convincingly explain what's going on in the experiment or
the analysis, keep reading or asking others until you can.
 
This package is meant to help you get started -- not to offload your
responsibility, or to substitute for your thinking.

### 2. Write everything down.

Meticulous notes are just as critical for computational work as they are for
experimental work. I keep records in three ways:
    
1. Scripts. These are text files containing the actual commands that perform
each analysis step. Scripts force you to keep track of your analyses.
 
2. README files. If there is anything important to know or remember about
what happened in a specific directory, I write it in a README text file
that I keep in that specific directory.
 
3. Electronic lab notebook. I keep a chronological lab notebook made of
text file per day (see [mdlabbook](http://github.com/fredpdavis/mdlabbook).
You should use whatever system you are comfortable with, but I highly
recommend you keep eletronic notes. When you want to remember what commands
you tried, or what analysis you were working on Monday two weeks ago, you
are unlikely to remember accurately unless you can look back at notes.

### 3. Trust no one, especially yourself

This is of course an exaggeration -- ultimately you are trusting multiple
layers of software that operate the sequencer through to the programs
that generate the figures you will intrepret. Be skeptical of all results.
Don't over-interpret or assume your results are correct: lots of things can
and do go wrong. Do you have positive and negative controls? If something
looks weird, it probably is.

### 4. Look inside.

Scripts are just text files -- view them in a text editor, or at the command
line using the `less FILENAME` command, to see the commands that are being
run.

### 5. Be tidy.

#### i. keep an organized directory structure.

I make a new folder for each project I work on, with the following structure:
    
- README -- file where I describe the project and broad goals
- data/ - directory storing raw data, each in its own folder -- eg, data/fastq/
- src/ - directory storing all scripts / code used in the project. I usually organize by language -- eg, src/R/
- run/ - directory where I run scripts, each kept in a dated directory -- eg, run/20181027.align_reads/
- results/ - store results from off-the-shelf programs -- eg, results/kallisto/
- analysis/ - where I do 'secondary' analysis, like making figures

This package comes with, or builds up, the directory structure above.

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

## RNA-seq

### Goal

We will focus on the most routine (and often most informative) analyses.

1. estimating transcript abundance (= expression level)

2. identifying differentially expressed genes between pairs of conditions

3. visualizing results using genome browsers, scatterplots and heatmaps.

We will not cover the many other kinds of analyses you can perform on RNA-seq
measurements, including identifying alternative splicing events, estimating
nascent transcription, or evaluating more complex experimental designs.

### Data sources

This package uses the following external data sources. __You don't need to download these__, as this package either comes with or will retrieve the appropriate files.

| data set                                                                                      | purpose                                       |
|-----------------------------------------------------------------------------------------------|-----------------------------------------------|
| [ENSEMBL GRCm38.94 release94](https://useast.ensembl.org/Mus_musculus/Info/Index)             | genome sequence, gene annotations             |
| [Life Tech: ERCC92](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip)        | synthetic spike-in sequences, concentratsion  |

### Software

This package uses the following underlying tools to analyze RNA-seq data.

| software                                                      | purpose                                                   |
|---------------------------------------------------------------|-----------------------------------------------------------|
| [kallisto 0.44.0](https://pachterlab.github.io/kallisto/)     | pseudoalign reads to transcriptome, estimates abundance   |
| [STAR 2.5.4a](https://github.com/alexdobin/STAR)              | align reads to genome                                     |
| [deeptools 3.1.2](https://deeptools.readthedocs.io/)          | make bigwig files for visualizing in IGV                  |
| [sleuth 0.29](https://pachterlab.github.io/sleuth/)           | R package for differential expression analysis            |
| [R 3.5.0](https://www.r-project.org)                          | creating figures, and running sleuth                      |

__You don't need to download these__, as they are already installed on the NIH HPC
and accessible through their very well maintained system of
['modules'](https://hpc.nih.gov/apps/modules.html)

These programs offer many different options -- the scripts in this package use a
particular incantation of these programs. As you get more comfortable with the
analysis, look inside the scripts and study the specific options used for each
program, and then read their manuals to understand what each does and what other
options you may want to try.


### 1. Setup your project directory

- Login to helix

    - on Mac: open the Terminal program, and type `ssh helix.nih.gov`, press enter, and enter your password when prompted

    - on Windows: download and run [PuTTY](https://www.putty.org) to ssh into helix.nih.gov, specifying your username and password

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
git clone https://github.com/fredpdavis/oshealab_genomic_analysis.git
```

- rename to your project name, cytokineX in this example:

```
mv oshealab_genomic_analysis cytokineX
```

- change directories to this project.

```
cd cytokineX
```

- List the files to see what you're getting in this package.

```
tree .
```

### 2. Prepare your input files

Now that you have a set of scripts and basic directory structure, you can start
using it to analyze your own data.

#### Get FASTQ sequence files

For testing, I provide four small fastq files-- you can delete if you want.

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

You can request your own sequence files from NIAMS core if your samples were
sequenced locally, or download them from
[NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) if you want to re-process
published data.

Store the sequence files in the `data/fastq/` directory in subdirectories
named by sequencing run identifiers (if sequenced locally), or by author if
the FASTQ comes from a published paper.


#### Edit the sample sheet `metadata/rnaseq_samples.txt`

- This is a tab-delimited text file. Make sure you use tab's, __NOT SPACES__

- This package includes an example file listing the test samples

```
shell> cat metadata/rnaseq_samples.txt

sampleID	sampleName	runID	fastqName	cytokineStim
r001	unstim_rep1	test	in1.fq.gz	nc
r002	unstim_rep3	test	in2.fq.gz	nc
r003	IFNg_72h_rep1	test	in3.fq.gz	IFNg_72h
r004	IFNg_72_rep2	test	in4.fq.gz	IFNg_72h
```

- Edit this file to describe your samples. nano is an easy-to-use text
editor: `nano metadata/rnaseq_samples.txt`

- requires 4 fields:
    1. sampleID
    2. sampleName
    3. runID
    4. fastqName

- add additional columns to describe sample properties like celltype,
  tissue, etc. this is useful for specifying the condition pairs you
  want to compare in the next file.

- expects to find fastq files in `data/fastq/<runID>/<fastqName>`


#### Edit the comparisons file, listing pairs of conditions to compare:

- This is a tab-delimited text file. Make sure you use tab's, __NOT SPACES__

- Each line represents a pairwise comparison

- Two columns, each specifying the samples to compare

- Samples can be specified by sampleName or by other features specified in
the sample sheet

```
shell> cat metadata/rnaseq_comparisons.txt

cytokineStim=nc	cytokineStim=IFNg_72h
```

- If you want to specify your samples using more than one feature, use
commas to express logical AND, and semicolon to express logical OR.
AND's will be interpreted first.

```
cytokineStim=nc,cellType=CD4	cytokineStim=nc,cellType=CD8
```

### 3. Retrieve and prepare genome and gene information

This step will retrieve the files necessary to process mouse data. Take a look
inside the script to get a sense of what data is being retrieved. If you ever
want to change genome version, gene annotations, species, etc. you will edit
this file.

- Make a new directory and copy the prepared_indices.slurm.csh script there.

```
mkdir -p run/20181102.prepare_inputs
cp src/slurm/prepare_indices.slurm.csh run/20181102.prepare_inputs
```

- Edit the script to specify the BASEDIR

```
set BASEDIR="/data/davisfp/projects/cytokineX"
```

- Submit the script to biowulf to download and process necessary files

```
ssh biowulf.nih.gov
cd /data/davisfp/projects/cytokineX/20181102.prepare_indices
mkdir slurm_out
sbatch prepare_indices.slurm.csh
```

- To check if the job is complete, use the squeue command

```
squeue -u davisfp
```

- Once the job is done, logoff biowulf

```
logout
```

You only need to run this step once (per project). You don't need to run this
step again if you just want to process additional samples.

### 4. Process your samples

- Make a new directory and copy the next script there.

```
mkdir -p run/20181102.process_samples
cp src/slurm/process_rnaseq_samples.slurm.csh run/20181102.process_samples
```

- Edit the script to specify the BASEDIR

```
set BASEDIR="/data/davisfp/projects/cytokineX"
```

- Edit the script to specify which tasks to process. By default, all samples
    listed in the metadata file will be processed. If you just want to process
    a subset of those samples, you can specify their defining features in the
    script
    
```
set SAMPLE_OPTION="-cytokine no -cellType CD4"
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

- Login to biowulf and submit the jobs to 

```
ssh biowulf.nih.gov
cd /data/davisfp/projects/cytokineX/20181102.process_samples
mkdir slurm_out
sbatch process_rnaseq_samples.slurm.csh
```

- Check job status with `squeue`; once the job is done, logoff biowulf

```
logout
```

### 5. Create figures and tables

The next steps of identifying differentially expresed genes and creating figures
is performed by an R script.

Unlike the shell scripts that we submitted to the cluster, I like to keep only
one R script, exactly where it sits in `src/R/basicRnaSeqAnalysis.R`

Edit this script (as you feel comfortable) to adapt it and change figures, etc.

We will run R on a cluster node -- instead of submitting a "batch" job with
sbatch, we will use sinteractive to request an interactive session. Thsi will
basically give you a command line, but one on a much bigger machine than your
desktop or laptop.

- login to biowulf and request an interactive session

```
ssh biowulf
cd data/projects/cytokineX
sinteractive --x11 --cpus-per-task=2 --mem=64g --time=24:00:00
```

- make an analysis directory

```
cd data/projects/cytokineX
mkdir -p analysis/20181102.makeFigures
cd analysis/20181102.makeFigures
```

- load the R 'module'

```
module load R
```

- start R and generate the figures

```
R> source("../../src/R/basicRnaSeqAnalysis.R")
R> dat <- loadData()
R> tx <- makeFigures(dat)
```

- This will create the figures in this directory

- To view the figures, secury copy (scp) the whole directory to your local desktop/laptop

```
scp -r data/projects/cytokineX/analysis/20181102.makeFigures .
```


## Design decisions

### Options

- genome version. default mm10
- annotation source/version. default ENSEMBL release 94
- ERCC spike-ins. default add to index. if not in sample, negligible contribution anyways and can always delete them.
- gene types. default only cDNA (protein-coding genes)
