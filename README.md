# oshealab_genomic_analysis -- scripts to perform routine genomics analyses

This package provides instructions and scripts to perform routine analysis of
genomics data (currently just RNA-seq) and encourage reproducibility.

Written for the HPC computing environment at NIH.

If you are new to computational work, please read the section on general advice.

## RNA-seq

### Goal

We will focus on the most routine (often most informative) analyses:
(1) estimating transcript abundance (= expression level)
(2) identifying differentially expressed genes between pairs of conditions
(3) visualizing results with scatterplots and heatmaps.

We will not cover the many other kinds of analyses you can perform on RNA-seq
measurements, including identifying alternative splicing events, estimating
nascent transcription, or evaluating more complex experimental designs.

### 1. Sample sheet

There are two tab-delimited text files that you should edit

__REMEMBER, these are tab-delimited, errant spaces will mess things up__

1. Sample description. 

### External software

software        purpose
----------      ---------------------------------------------------------------------
kallisto        pseudoaligns to transcriptome to estimate transcript abundance
STAR            aligns to genome (eg, for IGV visualization)
deeptools       make bigwig files
sleuth          R package for differential expression analysis
R               running sleuth, and creating figures


## General advice for computational work

1. Take responsibility.

    Neither the experiment nor the analysis is rocket science -- thousands of
    people have done these experiments before -- but, you should invest the time
    to fully understand the experiment and how its measurements are interpreted.
 
    These scripts are meant to help you get started -- not to offload your
    responsibility, or to substitute for your thinking.

2. Write everything down.

    It is critical to keep meticulous notes recording all details. Records are
    just as important for computational as it is for experimental work.

    I keep records in three ways:
    
    1. Scripts are text files containing commands that will be run. Organizing
       your analysis into scripts forces you to keep track of your analyses.
 
    2. README files -- If there is anything important to know or remember about
    what happened in a specific directory, I write it in a README text file
    that I keep in that specific directory.
 
    3. Electronic lab notebook. I keep a chronological lab notebook (one text
    file per day). You should use whatever system you are comfortable with,
    but I highly recommend you keep eletronic notes. When you want to remember
    a specific command you tried, or what analysis you were working on two
    weeks ago, you can flip back to your notes and figure it out.

3. Look inside.

    Scripts are simple text files -- open them in a text editor
    and look at the commands that are being run. It's not complicated.

4. Be skeptical.

    Don't over-interpret or take the results as immediately correct: lots of
    things can and do go wrong. Do you have positive and negative controls?
    If something looks weird, it probably is.

5. Be tidy.

    i. keep an organized directory structure.

    ii. don't use spaces or weird characters (asterisks, slashes, apostrophes,
       quotes, etc) in your sample names, condition names.

    iii. don't edit stuff unless you understand what you're doing.

    iv. Turn off your Mac's quote autocorrection. Apostrophes and quotes have
    specific meanings in scripts. By default, Mac's will often change these
    kinds of characters to curly quotes that will break the command -- this is
    a problem when eg, you try to test a command in a script by copying and
    pasting to a command line -- and then you find it fails.

## Organizing your directory

I make a new folder for each project I work on, with the following structure:

- README -- file where I describe the project and broad goals
- data/ - directory storing raw data, each in its own folder -- eg, data/fastq/
- src/ - directory storing all scripts / code used in the project. I usually organize by language -- eg, src/R/
- run/ - directory where I run scripts, each kept in a dated directory -- eg, run/20181027.align_reads/
- results/ - store results from off-the-shelf programs -- eg, results/kallisto/
- analysis/ - run scripts for 'secondary' analysis, like making plots

