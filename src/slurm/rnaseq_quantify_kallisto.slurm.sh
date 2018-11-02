#!/bin/csh -v
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-core=1
#SBATCH --mem=10g
#SBATCH -o slurm_out/rnaseq_quantify_kallisto.%A.%a.out
#SBATCH -e slurm_out/rnaseq_quantify_kallisto.%A.%a.err
#SBATCH --time=10:00:00
#SBATCH --gres=lscratch:200
#SBATCH --array=1-33

set NUMCPU=8

set BASEDIR="/data/davisfp/projects/cytokineX"
set SAMPLEINFO_FN="$BASEDIR/metadata/rnaseq_samples.txt"
set SAMPLE_OPTION=""

set SCRATCHDIR="/lscratch/$SLURM_JOBID"

# Module load
module load kallisto/0.42.4
module load STAR/2.5.4a

set SRCDIR="$BASEDIR/src"
set names=( `perl $SRCDIR/perl/txt2tasklist.pl $SAMPLEINFO_FN sampleName $SAMPLE_OPTION` )
set runids=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN runID $SAMPLE_OPTION` )
set fastqs=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN fastqFn $SAMPLE_OPTION` )

set SAMPLENAME=$names[$SLURM_ARRAY_TASK_ID]
set FASTQ=$fastqs[$SLURM_ARRAY_TASK_ID]
set RUNID=$runids[$SLURM_ARRAY_TASK_ID]

### SETUP NECESSARY DIRECTORIES
set FASTQ_FN="$BASEDIR/data/fastq/$RUNID/$FASTQ"

set KALLISTO_OUT_BASEDIR="$BASEDIR/results/RNAseq/kallisto/$RUNID"
set KALLISTO_OUTDIR="$KALLISTO_OUT_BASEDIR/$SAMPLENAME"
set KALLISTO_INDEX="$BASEDIR/data/kallisto_files/GRCm38.82.FPtags.ERCC.kallisto_index"

### SETUP NECESSARY FILES
##set KALLISTO_BIN="kallisto"
set KALLISTO_BIN="/data/davisfp/software/kallisto/kallisto_linux-v0.42.4/kallisto_linux-v0.42.4/kallisto"

foreach T_OUTDIR ( $KALLISTO_OUTDIR )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end


### SETUP PROGRAM OPTIONS
set KALLISTO_OPTIONS="quant -t $NUMCPU -i $KALLISTO_INDEX --single -l 200 -s 30 -b 50 -o $KALLISTO_OUTDIR $FASTQ_FN"

### START ACTUALLY RUNNING
set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#SLURM run started on $curhost at $curtime"
echo "#scratchdir = $SCRATCHDIR"

echo "# $KALLISTO_BIN $KALLISTO_OPTIONS ($curtime)"
$KALLISTO_BIN $KALLISTO_OPTIONS

set curtime=`date`
echo "#SLURM run finished on $curhost at $curtime"
