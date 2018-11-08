#!/bin/csh -v
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-core=1
#SBATCH --mem=100g
#SBATCH -o slurm_out/process_rnaseq_samples.%A.%a.out
#SBATCH -e slurm_out/process_rnaseq_samples.%A.%a.err
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:200
#SBATCH --array 1-4


set NUMCPU=8
set scratchdir="/lscratch/$SLURM_JOBID"

set BASEDIR="/data/davisfp/projects/cytokineX"
set SAMPLEINFO_FN="$BASEDIR/metadata/rnaseq_samples.txt"
set SAMPLE_OPTION=""

set sampleids=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN sampleID $SAMPLE_OPTION` )
set samplenames=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN sampleName $SAMPLE_OPTION` )
set runids=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN runID $SAMPLE_OPTION` )
set fastqs=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN fastqName $SAMPLE_OPTION` )

set SAMPLEID=$sampleids[$SLURM_ARRAY_TASK_ID]
set SAMPLENAME=$samplenames[$SLURM_ARRAY_TASK_ID]
set RUNID=$runids[$SLURM_ARRAY_TASK_ID]
set FASTQ=$fastqs[$SLURM_ARRAY_TASK_ID]

## LOAD PROGRAMS
module load kallisto/0.44.0
module load STAR/2.5.4a
module load picard/2.17.11
module load java/1.8.0_92
module load R/3.5.0
module load deeptools/3.1.2

set PICARDSORT_BIN="java -jar $PICARDJAR SortSam"
set PICARDINDEX_BIN="java -jar $PICARDJAR BuildBamIndex"
set PICARDSTATS_BIN="java -jar $PICARDJAR CollectRnaSeqMetrics"

## SPECIFY INPUT FILES
set ENSEMBL_GENOMEVER="GRCm38"
set ENSEMBL_RELEASE="94"
set EGR="${ENSEMBL_GENOMEVER}.${ENSEMBL_RELEASE}"
set KALLISTO_INDEX="$BASEDIR/data/kallisto_files.$EGR/$EGR.kalisto_index"
set STARDIR="$BASEDIR/data/star_files.$EGR"
set STAR_INDEX="$STARDIR/$EGR.star_index"
set UNCOMPR_GTF="$STARDIR/Mus_musculus.$EGR.ERCC92.gtf"
set PICARD_REF_FLAT_FN="$BASEDIR/data/picard_files.$EGR/refFlat.$EGR.ERCC92.txt.gz"

set FASTQ_FN="$BASEDIR/data/fastq/$RUNID/$FASTQ"

## OUTPUT DIRECTORIES
set KALLISTO_OUTDIR="$BASEDIR/results/RNAseq/kallisto.$EGR/$SAMPLEID"
set STAR_OUTDIR="$BASEDIR/results/RNAseq/star.$EGR/$SAMPLEID"
set PICARDSTATS_OUTDIR="$BASEDIR/results/RNAseq/picard_stats.star.$EGR/$SAMPLEID"

foreach T_OUTDIR ( $KALLISTO_OUTDIR $STAR_OUTDIR $PICARDSTATS_OUTDIR )
   if (! -e $T_OUTDIR) then
      echo "mkdir -p $T_OUTDIR"
      mkdir -p $T_OUTDIR
   endif
end


### OUTPUT FILES
set STAR_OUT_BAM_FN="${STAR_OUTDIR}/Aligned.out.bam"
set STAR_OUT_SORTED_BAM_FN="${STAR_OUTDIR}/Aligned.out.sorted.bam"
set STAR_OUT_SORTED_BAI_FN="${STAR_OUTDIR}/Aligned.out.sorted.bai"
set STAR_OUT_BW_FN="${STAR_OUTDIR}/$SAMPLEID.star_deeptool_rpkm.bw"
set STAR_NAMED_OUT_BW_FN="${STAR_OUTDIR}/${SAMPLENAME}.star_deeptool_rpkm.bw"

set PICARDSTATS_TXT_FN="${PICARDSTATS_OUTDIR}/${SAMPLEID}.picard_rnaseq_report.txt"
set PICARDSTATS_PDF_FN="${PICARDSTATS_OUTDIR}/${SAMPLEID}.picard_rnaseq_report.pdf"

### PROGRAM OPTIONS
set KALLISTO_OPTIONS="quant -t $NUMCPU -i $KALLISTO_INDEX --single -l 200 -s 30 -b 50 -o $KALLISTO_OUTDIR $FASTQ_FN"
set STAR_OPTIONS="--genomeDir $STAR_INDEX --readFilesIn $FASTQ_FN --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterType BySJout --runThreadN $NUMCPU --outFileNamePrefix $STAR_OUTDIR/ --readFilesCommand zcat"

set PICARDSORT_OPTIONS="I=$STAR_OUT_BAM_FN O=$STAR_OUT_SORTED_BAM_FN SO=coordinate"
set PICARDINDEX_OPTIONS="I=$STAR_OUT_SORTED_BAM_FN O=$STAR_OUT_SORTED_BAI_FN"
set PICARDSTATS_OPTIONS="REF_FLAT=$PICARD_REF_FLAT_FN STRAND_SPECIFICITY=NONE INPUT=$STAR_OUT_SORTED_BAM_FN CHART_OUTPUT=$PICARDSTATS_PDF_FN OUTPUT=$PICARDSTATS_TXT_FN"

set BAMCOVERAGE_OPTIONS="-b $STAR_OUT_SORTED_BAM_FN -p $NUMCPU -o $STAR_OUT_BW_FN --normalizeUsing RPKM"


### START ACTUALLY RUNNING

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`

echo "# slurm run started on $curhost at $curtime"
echo "# Processing sampleID $SAMPLEID ($SAMPLENAME)" ;

set echo

set curtime=`date`
echo "STEP 1. KALLISTO PSEUDOALIGNMENT ($curtime)"
kallisto $KALLISTO_OPTIONS

set curtime=`date`
echo "STEP 2. STAR ALIGNMENT, PICARD SORT&INDEX ($curtime)"
STAR $STAR_OPTIONS
$PICARDSORT_BIN $PICARDSORT_OPTIONS
$PICARDINDEX_BIN $PICARDINDEX_OPTIONS

set curtime=`date`
echo "STEP 3. DEEPTOOLS GENOME TRACKS ($curtime)"
bamCoverage $BAMCOVERAGE_OPTIONS
ln -s $STAR_OUT_BW_FN $STAR_NAMED_OUT_BW_FN

set curtime=`date`
echo "STEP 4. PICARD QC ($curtime)"
$PICARDSTATS_BIN $PICARDSTATS_OPTIONS

set curtime=`date`
echo "# slurm run finished on $curhost at $curtime"
