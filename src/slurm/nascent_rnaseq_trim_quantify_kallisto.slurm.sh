#!/bin/csh -v
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-core=1
#SBATCH --mem=10g
#SBATCH -o slurm_out/rnaseq_quantify_kallisto.%A.%a.out
#SBATCH -e slurm_out/rnaseq_quantify_kallisto.%A.%a.err
#SBATCH --time=10:00:00
#SBATCH --gres=lscratch:200
#SBATCH --array=1-9

set NUMCPU=2
set SCRATCHDIR="/lscratch/$SLURM_JOBID"
set BASEDIR="/gpfs/gsfs4/users/davisfp/projects/exfoxp3"

set SAMPLEINFO_FN="$BASEDIR/data/exfoxp3_rna_samples.txt"
set SAMPLE_OPTION="dataType2 nascent"

# Module load
module load kallisto/0.43.0
module load trimmomatic/0.36

set SRCDIR="$BASEDIR/src"
set names=( `perl $SRCDIR/perl/txt2tasklist.pl $SAMPLEINFO_FN sampleName $SAMPLE_OPTION` )
set runids=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN runID $SAMPLE_OPTION` )
set fastqs=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN fastqFn $SAMPLE_OPTION` )

set SAMPLENAME=$names[$SLURM_ARRAY_TASK_ID]
set FASTQ=$fastqs[$SLURM_ARRAY_TASK_ID]
set RUNID=$runids[$SLURM_ARRAY_TASK_ID]

### SETUP NECESSARY DIRECTORIES
set FASTQ_FN="$BASEDIR/data/fastq/$RUNID/$FASTQ"

set KALLISTO_OUT_BASEDIR="$BASEDIR/results/RNAseq/kallisto_stranded/$RUNID"
set KALLISTO_OUTDIR="$KALLISTO_OUT_BASEDIR/$SAMPLENAME"
set KALLISTO_INDEX="$BASEDIR/data/kallisto_files/GRCm38.82.FPtags.ERCC.kallisto_index"

set TRIM_FASTQ_DIR="$BASEDIR/results/RNAseq/trimmed_fastq/$RUNID"
set TRIM_FASTQ_FN="$TRIM_FASTQ_DIR/trim_$FASTQ"

### SETUP NECESSARY FILES
set TRIMMOMATIC_BIN="java -jar $TRIMMOJAR"
set KALLISTO_BIN="kallisto"

foreach T_OUTDIR ( $TRIM_FASTQ_DIR $KALLISTO_OUTDIR )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end


### SETUP PROGRAM OPTIONS
set TRIMMOMATIC_OPTIONS="SE $FASTQ_FN $TRIM_FASTQ_FN HEADCROP:3"
set KALLISTO_OPTIONS="quant -t $NUMCPU -i $KALLISTO_INDEX --single -l 200 -s 30 --fr-stranded -b 50 -o $KALLISTO_OUTDIR $TRIM_FASTQ_FN"

### START ACTUALLY RUNNING
set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#SLURM run started on $curhost at $curtime"
echo "#scratchdir = $SCRATCHDIR"

if (! -e $TRIM_FASTQ_FN) then
   echo "# $TRIMMOMATIC_BIN $TRIMMOMATIC_OPTIONS ($curtime)"
   $TRIMMOMATIC_BIN $TRIMMOMATIC_OPTIONS
endif

echo "# $KALLISTO_BIN $KALLISTO_OPTIONS ($curtime)"
$KALLISTO_BIN $KALLISTO_OPTIONS

set curtime=`date`
echo "#SLURM run finished on $curhost at $curtime"
