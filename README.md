# oshealab_genomic_analysis -- scripts to perform routine genomics analyses

This package will help you get started analyzing genomics data (currently
only RNA-seq) in a scripted and reproducible way.

Written for the HPC computing environment at NIH.

## General advice for computational work

1. Take responsibility.

    Neither the experiment nor the analysis is complicated -- thousands of
    people have done RNA-seq before -- but, you still have to invest time to
    fully understand the experiment and how its measurements are interpreted.
    If you can't convincingly explain what's going on, keep reading or
    asking others until you can.
 
    This package is meant to help you get started -- not to offload your
    responsibility, or to substitute for your thinking.

2. Write everything down.

    Meticulous notes are just as critical for computational work as it is for
    experimental work.  I keep records in three ways:
    
    1. Scripts are text files containing commands that will be run. Organizing
       your analysis into scripts forces you to keep track of your analyses.
 
    2. README files -- If there is anything important to know or remember about
    what happened in a specific directory, I write it in a README text file
    that I keep in that specific directory.
 
    3. Electronic lab notebook. I keep a chronological lab notebook through a
    series of text files -- one per day. You should use whatever system you are
    comfortable with, but I highly recommend you keep eletronic notes. When you
    want to remember a specific command you tried, or what analysis you were
    working on two weeks ago, you are unlikely to remember it accurately unless
    you can look back at to your notes.

3. Look inside.

    Scripts are just text files -- view them in a text editor, or at the command
    line using the `less FILENAME` command, to see the commands that are being
    run.

4. Be skeptical.

    Don't over-interpret or take your results as immediately correct: lots of
    things can and do go wrong. Do you have positive and negative controls?
    If something looks weird, it probably is.

5. Be tidy.

    i. keep an organized directory structure.

        I make a new folder for each project I work on, with the following structure:
        
        - README -- file where I describe the project and broad goals
        - data/ - directory storing raw data, each in its own folder -- eg, data/fastq/
        - src/ - directory storing all scripts / code used in the project. I usually organize by language -- eg, src/R/
        - run/ - directory where I run scripts, each kept in a dated directory -- eg, run/20181027.align_reads/
        - results/ - store results from off-the-shelf programs -- eg, results/kallisto/
        - analysis/ - where I do 'secondary' analysis, like making figures

    ii. don't use spaces or weird characters (asterisks, slashes, apostrophes,
       quotes, etc) in your sample names, condition names. Use underscores and
       periods if you need separators.

    iii. don't edit stuff unless you understand what you're doing.

    iv. Turn off Mac autocorrect for quotes. Apostrophes and quotes have
    specific meanings in scripts. By default, Mac's will often change these
    kinds of characters to curly quotes that will break the command -- this is
    a problem when eg, you try to copy and paste a command into the terminal
    to test -- because the command line (or shell) doesn't like fancy quotes.

## RNA-seq

### Goal

We will focus on the most routine (often most informative) analyses:
(1) estimating transcript abundance (= expression level)
(2) identifying differentially expressed genes between pairs of conditions
(3) visualizing results with scatterplots and heatmaps.

We will not cover the many other kinds of analyses you can perform on RNA-seq
measurements, including identifying alternative splicing events, estimating
nascent transcription, or evaluating more complex experimental designs.


### 0. Setup your project directory

You will setup this analysis on the helix machine.

- Login to helix, download this package, and rename the directory to your
  project name  (cytokineX in this example).

```sh
yourlaptop> ssh helix.nih.gov
helix> cd /data/davisfp/projects
helix> git clone https://github.com/fredpdavis/oshealab_genomic_analysis.git
helix> mv oshealab_genomic_analysis cytokineX
```

### 1. Prepare your input files

1. FASTQ sequence files -- `data/fastq`

    You can request these from the NIAMS core (or download them from
    [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) if you want to re-process
    published data.

    Store these in the `data/fastq/` directory in subdirectories named by
    either flowcell, if they were locally sequenced, or by author if the
    FASTQ comes from a published paper.

    These files often have a .fastq or .fq suffix, and are typically compressed
    by either gzip or bunzip, resulting in a further '.gz' or '.bz2' suffix.

    I provide four short fastq files for testing -- you can delete these if you want.

    ```
    ls data/fastq/test

    in1.fq.gz
    in2.fq.gz
    ```

2. Sample Sheet: `metadata/rnaseq_samples.txt`

    - tab-delimited text file

    - Edit this file in the metadata directory to describe your samples
        - nano is an easy-to-use text editor that you can invoke at the
            command line: `nano metadata/rnaseq_samples.txt`

    - requires 3 fields:
        1. flowcell
        2. fastqName
        3. sampleName

    - add additional columns to describe sample properties like celltype,
      tissue, etc. this is useful for specifying the condition pairs you
      want to compare in the next file.

    - expects to find fastq files in BASEDIR/data/fastq/<flowcell>/<fastqName>

    ```
    cat metadata/rnaseq_samples.txt
    flowcell	fastqName	sampleName	cytokineStim
    test	in1.fq.gz	unstim_rep1	no
    test	in2.fq.gz	unstim_rep2	no
    test	in3.fq.gz	cytokineX_stim_1hr_rep1	cytokineX_1hr
    test	in4.fq.gz	cytokineX_stim_1hr_rep2	cytokineX_1hr
    ```

3. List of comparisons -- list of condition pairs to compare

    - tab-delimited text file
    - Each line represents a pairwise comparison
    - Two columns, each specifying the samples to compare
    - Samples can be specified by sampleName or by other features specified in
      the sample sheet

    ```
    > cat metadata/rnaseq_comparisons.txt
    cytokineStim=no	cytokineStim=cytokineX_1hr
    ```

    - If you want to specify your samples using more than one feature, use
    commas to express logical AND, and semicolon to express logical OR.
    AND's will be interpreted first.

    ```
    cytokineStim=no,cellType=CD4	cytokineStim=no,cellType=CD8
    ```

__REMEMBER, these files are tab-delimited, errant spaces will mess things up__

### 2. Retrieve and prepare genome/gene information

This step will retrieve the files necessary to process mouse data. Take a look
inside the script to get a sense of what data is being retrieved. If you ever
want to change genome version, gene annotations, species, etc. you will edit
this file.

- Make a new directory and copy the script there.

    ```sh
    mkdir -p run/20181102.prepare_inputs
    cp src/slurm/prepare_input_files.slurm.csh run/20181102.prepare_inputs
    ```

- Edit the script to specify the BASEDIR

    ```sh
    set BASEDIR="/data/davisfp/projects/cytokineX"
    ```

- Submit the jobs to biowulf to download and process necessary files

    ```sh
    ssh biowulf.nih.gov
    cd /data/davisfp/projects/cytokineX/20181102.prepare_inputs
    mkdir slurm_out
    sbatch prepare_input_files.slurm.csh
    ```

- To check if the job is complete, use the squeue command

- Once the job is done, logoff biowulf

    ```sh
    logout
    ```

You only need to run this step once (per project). You don't need to run this
step again if you just want to process additional samples.

### 2. Run the primary processing.

- Make a new directory and copy the next script there.

    ```sh
    mkdir -p run/20181102.align_samples
    cp src/slurm/align_samples.slurm.csh run/20181102.prepare_inputs
    ```

- Edit the script to specify the BASEDIR

    ```sh
    set BASEDIR="/data/davisfp/projects/cytokineX"
    ```

- Edit the script to specify which tasks to process. By default, all samples
    listed in the metadata file will be processed. If you just want to process
    a subset of those samples, you can specify their defining features in the
    script
    
    ```sh
    set SAMPLE_OPTION="-cytokine no -cellType CD4"
    ```

- Edit the script to tell the cluster how many jobs you will run. For example,
  we'd like to process the 4 samples in our test set:

    ```sh
    #SBATCH --array=1-4
    ```

- It makes sense to run this script on a single job to begin with, just to
  make sure that everything runs properly:

    ```sh
    #SBATCH --array=1-1
    ```

- Login to biowulf and submit the jobs to 

    ```sh
    ssh biowulf.nih.gov
    cd /data/davisfp/projects/cytokineX/20181102.align_samples
    mkdir slurm_out
    sbatch align_samples.slurm.csh
    ```

- Check job status with `squeue`; once the job is done, logoff biowulf

    ```sh
    logout
    ```

### 2. Run the secondary processing

The next steps of identifying differentially expresed genes and creating figures
is performed by an R script.

Unlike the shell scripts that we submitted to the cluster, I like to keep only
one R script, exactly where it sits in `src/R/basicRnaAnalysis.R`

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
    R
    R> source("../../src/R/basicRnaAnalysis.R")
    R> dat <- loadData()
    R> tx <- makeFigures(dat)
    ```

- This will create the figures in this directory

- To view the figures, copy the whole directory to your local desktop/laptop

    ```
    scp -r data/projects/cytokineX/analysis/20181102.makeFigures .
    ```

### External software

| software   |    purpose                                                      |
|------------|-----------------------------------------------------------------|
| kallisto   |    pseudoaligns reads to transcriptome, estimates abundance     |
| STAR       |    aligns reads to genome                                       |
| deeptools  |    make bigwig files for visualizing in IGV                     |
| sleuth     |    R package for differential expression analysis               |
| R          |    running sleuth, and creating figures                         |
