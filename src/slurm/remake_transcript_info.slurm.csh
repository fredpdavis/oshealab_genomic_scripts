#!/bin/csh -v

set BASEDIR="/data/davisfp/projects/cytokineX"
set NUMCPU=8

## SPECIFY INFORMATION ABOUT ASSEMBLy 
set ENSEMBL_GENOMEVER="GRCm38"
set ENSEMBL_RELEASE="94"
set EGR="${ENSEMBL_GENOMEVER}.${ENSEMBL_RELEASE}"
set GENOMEDIR="$BASEDIR/data/ENSEMBL.$EGR"
set GENOMEFA="$GENOMEDIR/genome.$EGR.ERCC92.fa"
set TRANSCRIPTFA="$GENOMEDIR/Mus_musculus.${ENSEMBL_GENOMEVER}.cdna.all.fa.gz"
set ERCCGTF="$BASEDIR/data/ercc/ERCC92.gtf"

set KALLISTODIR="$BASEDIR/data/kallisto_files.$EGR"
set TXINFO_FN="$KALLISTODIR/transcript_info.$EGR.txt"

set ORIG_GTF_FN="$GENOMEDIR/Mus_musculus.$EGR.gtf.gz"
set WITHERCC_GTF_FN="Mus_musculus.$EGR.ERCC92.gtf"

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`

echo "# slurm run started on $curhost at $curtime"

cd $KALLISTODIR
zcat $ORIG_GTF_FN > $WITHERCC_GTF_FN
cat $ERCCGTF >> $WITHERCC_GTF_FN
perl $BASEDIR/src/perl/GTF2transcript_info.pl < $WITHERCC_GTF_FN > $TXINFO_FN

cd $curdir

set curtime=`date`
echo "# slurm run finished on $curhost at $curtime"
