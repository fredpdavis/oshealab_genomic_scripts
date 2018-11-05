#!/bin/csh -v
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-core=1
#SBATCH --mem=100g
#SBATCH -o slurm_out/rnaseq_align_STAR.%A.%a.out
#SBATCH -e slurm_out/rnaseq_align_STAR.%A.%a.err
#SBATCH --time=6:00:00
#SBATCH --gres=lscratch:200
#SBATCH --array=1-9

set NUMCPU=15
set scratchdir="/lscratch/$SLURM_JOBID"
set BASEDIR="/data/davisfp/projects/exfoxp3"

set SAMPLEINFO_FN="$BASEDIR/metadata/exfoxp3_rna_samples.txt"
set SAMPLE_OPTION="dataType2 nascent"

# Module load
module load STAR/2.5.4a
module load samtools/0.1.19
module load bedtools/2.19.1
module load ucsc/365
module load deeptools

set SRCDIR="$BASEDIR/src"
set names=(  `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN sampleName $SAMPLE_OPTION` )
set runids=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN runID $SAMPLE_OPTION` )
set fastqs=( `perl $BASEDIR/src/perl/txt2tasklist.pl $SAMPLEINFO_FN fastqFn $SAMPLE_OPTION` )

set SAMPLENAME=$names[$SLURM_ARRAY_TASK_ID]
set FASTQ=$fastqs[$SLURM_ARRAY_TASK_ID]
set RUNID=$runids[$SLURM_ARRAY_TASK_ID]

### SETUP NECESSARY DIRECTORIES
set FASTQ_FN="$BASEDIR/data/fastq/$RUNID/$FASTQ"

set STAR_OUT_BASEDIR="$BASEDIR/results/RNAseq/STAR/$RUNID"
set STAR_OUTDIR="$STAR_OUT_BASEDIR/$SAMPLENAME"
set STAR_INDEX="$BASEDIR/data/star_files/GRCm38.82.ERCC.star_index"

set STAR_OUT_SAM_FN="${STAR_OUTDIR}/Aligned.out.sam"
set STAR_OUT_BAM_PREFIX="${STAR_OUTDIR}/${SAMPLENAME}.sorted"
set STAR_OUT_BAM_FN="${STAR_OUT_BAM_PREFIX}.bam"
set STAR_OUT_BAI_FN="${STAR_OUT_BAM_PREFIX}.bam.bai"
set STAR_OUT_BED_FN="${STAR_OUTDIR}/${SAMPLENAME}.star_scaled10M_nonRRNA.bed"
set STAR_OUT_BW_FN="${STAR_OUTDIR}/${SAMPLENAME}.star_scaled10M_nonRRNA.bw"
set STAR_FWD_BW_FN="${STAR_OUTDIR}/${SAMPLENAME}.star_fwd.bw"
set STAR_REV_BW_FN="${STAR_OUTDIR}/${SAMPLENAME}.star_rev.bw"

### SETUP NECESSARY FILES
set STAR_BIN="STAR"
set SAMTOOLS_BIN="samtools"
set CONVERT_ENSEMBL_TO_UCSC_BED_BIN="perl $SRCDIR/perl/convert_ensembl_to_ucsc_chrom.pl"
set COUNTFASTQ_BIN="$SRCDIR/perl/count_fastq_reads.pl"
set SCALEBED_BIN="$SRCDIR/perl/scale_bed_rpm.pl"
set WIGTOBIGWIG_BIN="wigToBigWig"
set TEMP_CHROMSIZES_FN="${STAR_OUTDIR}/bam_chromsizes.txt"
set TEMP_CHROMSIZES_UCSC_FN="${STAR_OUTDIR}/bam_chromsizes_ucsc.txt"

set ALICOUNT_FN="$SAMPLENAME.alignment_stats.txt"
set RRNA_BED_FN="$BASEDIR/data/star_files/rRNA_regular_chromosomes.bed"
set NONRRNA_BED_FN="$BASEDIR/data/star_files/regular_chromosomes.nonribosomal_regions.bed"


foreach T_OUTDIR ( $STAR_OUTDIR $scratchdir )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end


### SETUP PROGRAM OPTIONS
set STAR_OPTIONS="--genomeDir $STAR_INDEX --readFilesIn $FASTQ_FN --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterType BySJout --runThreadN $NUMCPU --outFileNamePrefix $STAR_OUTDIR/ --readFilesCommand zcat"

### START ACTUALLY RUNNING
set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#sgejob run started on $curhost at $curtime"
echo "# Procesing: $SAMPLENAME ($RUNID)" ;

if (! -e "${STAR_OUTDIR}/Aligned.out.sam" ) then
   set curtime=`date`
   echo "# Running STAR: $STAR_BIN $STAR_OPTIONS ( $curtime )"
   $STAR_BIN $STAR_OPTIONS
endif

if (! -e $STAR_OUT_BAM_FN ) then
   set curtime=`date`
   echo "# Converting to sorted BAM file: $SAMTOOLS_BIN view -@ $NUMCPU -b -S $STAR_OUT_SAM_FN | $SAMTOOLS_BIN sort -@ $NUMCPU -m 5G - $STAR_OUT_BAM_PREFIX ($curtime)"
   $SAMTOOLS_BIN view -@ $NUMCPU -b -S $STAR_OUT_SAM_FN | $SAMTOOLS_BIN sort -@ $NUMCPU -m 5G - $STAR_OUT_BAM_PREFIX
   unlink $STAR_OUT_SAM_FN
endif

if (! -e "${STAR_OUTDIR}/$ALICOUNT_FN") then
   set curtime=`date`
   echo "# Counting alignment statistics ($curtime)"
   set NUMREADS=`perl $COUNTFASTQ_BIN $FASTQ_FN`
   echo "# Number of reads: $NUMREADS" > ${STAR_OUTDIR}/$ALICOUNT_FN

   set NUMALN_READS=`$SAMTOOLS_BIN view $STAR_OUT_BAM_FN | cut -f1 | sort -u -T${scratchdir} | wc -l`
   echo "# Number of aligned reads: $NUMALN_READS" >> ${STAR_OUTDIR}/$ALICOUNT_FN

   set NUMALN_RRNA_READS=`intersectBed -abam $STAR_OUT_BAM_FN -b $RRNA_BED_FN | $SAMTOOLS_BIN view - | cut -f1 | sort -u -T${scratchdir} | wc -l`
   echo "# Number of reads aligned to rRNA loci: $NUMALN_RRNA_READS" >> ${STAR_OUTDIR}/$ALICOUNT_FN

   set NUMALNS=`$SAMTOOLS_BIN flagstat $STAR_OUT_BAM_FN | grep -m 1 mapped | awk '{print $1}'`
   echo "# Number of alignments: $NUMALNS" >> ${STAR_OUTDIR}/$ALICOUNT_FN

   set NUMALNS_NONRRNA=`intersectBed -abam $STAR_OUT_BAM_FN -b $NONRRNA_BED_FN | $SAMTOOLS_BIN view -c -`
   echo "# Number of non-ribosomal, non-ERCC, non-INTACT alignments: $NUMALNS_NONRRNA" >> ${STAR_OUTDIR}/$ALICOUNT_FN
endif

$SAMTOOLS_BIN view -H $STAR_OUT_BAM_FN | grep "@SQ" | sed -e "s/SN://" -e "s/LN://" | cut -f 2,3 > $TEMP_CHROMSIZES_FN
$CONVERT_ENSEMBL_TO_UCSC_BED_BIN < $TEMP_CHROMSIZES_FN > $TEMP_CHROMSIZES_UCSC_FN

if (! -e $STAR_OUT_BW_FN ) then

   if (! -e $STAR_OUT_BED_FN ) then
      set curtime=`date`
      echo "# Creating normalized bedgraph of STAR results ( $curtime )"

      set NUMALNS=`$SAMTOOLS_BIN flagstat $STAR_OUT_BAM_FN | grep -m 1 mapped | awk '{print $1}'`
      genomeCoverageBed -split -bg -ibam $STAR_OUT_BAM_FN -g ${TEMP_CHROMSIZES_FN} | perl $SCALEBED_BIN $NUMALNS_NONRRNA 10000000 > $STAR_OUT_BED_FN
   endif

   set curtime=`date`
   echo "# Creating normalized bigWig of STAR results ( $curtime )"
   cat $STAR_OUT_BED_FN | $CONVERT_ENSEMBL_TO_UCSC_BED_BIN | $WIGTOBIGWIG_BIN -clip stdin $TEMP_CHROMSIZES_UCSC_FN $STAR_OUT_BW_FN 
   unlink $STAR_OUT_BED_FN

endif

# make stranded bigwig library
if (! -e $STAR_FWD_BW_FN ) then
   set curtime=`date`
   echo "# Creating stranded normalized bigWig of STAR results ( $curtime )"

   if (! -e $STAR_OUT_BAI_FN) then
      $SAMTOOLS_BIN index $STAR_OUT_BAM_FN
   endif

   bamCoverage --normalizeUsing RPKM -b $STAR_OUT_BAM_FN -o $STAR_FWD_BW_FN --samFlagExclude 16 -p $NUMCPU
   bamCoverage --normalizeUsing RPKM -b $STAR_OUT_BAM_FN -o $STAR_REV_BW_FN --samFlagInclude 16 -p $NUMCPU
endif

##unlink $STAR_OUT_BAM_FN #don't need the BAM file, delete it.
unlink $TEMP_CHROMSIZES_FN
unlink $TEMP_CHROMSIZES_UCSC_FN

set curtime=`date`
echo "#sgejob run finished on $curhost at $curtime"
