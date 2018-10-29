# oshealab_genomic_analysis -- scripts to perform routine genomics analyses

This package will help you get started analyzing genomics data (currently
only RNA-seq) in a reproducible way.

Written for the HPC computing environment at NIH.

## General advice for computational work

1. Take responsibility.

    Neither the experiment nor the analysis is complicated -- thousands of
    people have done this before -- but, you have to invest time to fully
    understand the experiment and how its measurements are interpreted.
    If you can't convincingly explain what's going on, keep reading or
    asking others until you can.
 
    This package is meant to help you get started -- not to offload your
    responsibility, or to substitute for your thinking.

2. Write everything down.

    Meticulous notes are just as critical for computational work as it is for
    experimental work.

    I keep records in three ways:
    
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

    Scripts are just text files -- open them in a text editor and look at the
    commands that are being run. It's not complicated.

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
        - analysis/ - run scripts for 'secondary' analysis, like making plots

    ii. don't use spaces or weird characters (asterisks, slashes, apostrophes,
       quotes, etc) in your sample names, condition names.

    iii. don't edit stuff unless you understand what you're doing.

    iv. Turn off Mac autocorrect for quotes. Apostrophes and quotes have
    specific meanings in scripts. By default, Mac's will often change these
    kinds of characters to curly quotes that will break the command -- this is
    a problem when eg, you try to copy and paste a command into the terminal
    to test -- because, it will fail.

## RNA-seq

### Goal

We will focus on the most routine (often most informative) analyses:
(1) estimating transcript abundance (= expression level)
(2) identifying differentially expressed genes between pairs of conditions
(3) visualizing results with scatterplots and heatmaps.

We will not cover the many other kinds of analyses you can perform on RNA-seq
measurements, including identifying alternative splicing events, estimating
nascent transcription, or evaluating more complex experimental designs.

### 1. Prepare your input files

This pipeline expects two text files as input.

1. Sample Sheet -- list of individual samples

    - First line is the header, which includes column names
    - requires 3 fields:
        1. flowcell
        2. fastqName
        3. sampleName

    - add additional columns to describe sample properties like celltype,
      tissue, etc. this is useful for specifying the condition pairs you
      want to compare in the next file.

    ```sh
    flowcell	fastqName	sampleName
    ```

    - expects to find the fastq file in BASEDIR/data/<flowcell>/<fastqName>

2. List of comparisons -- list of condition pairs to compare

    - Each line represents a pairwise comparison
    - Two columns, each specifying the samples to compare
    - Samples can be specified by sampleName or by other features specified in
      the sample sheet

    ```sh
    celltype=Foxp3,tissue=lung	celltype=Foxp3,tissue=spleen
    ```

- commas are interpreted as logical AND's; semicolons are interpreted as OR's
- logical order: AND's will be interpreted first

__REMEMBER, these files are tab-delimited, errant spaces will mess things up__

### 2. Prepare software-specific input files

You only need to run this step once. If you generate a new library, you don't
need to run this step again.

### 2. Run the primary processing.

### 3. Run the secondary processing: figures, tables.

### External software

| software   |    purpose                                                      |
|------------|-----------------------------------------------------------------|
| kallisto   |    pseudoaligns reads to transcriptome, estimates abundance     |
| STAR       |    aligns reads to genome                                       |
| deeptools  |    make bigwig files for visualizing in IGV                     |
| sleuth     |    R package for differential expression analysis               |
| R          |    running sleuth, and creating figures                         |
