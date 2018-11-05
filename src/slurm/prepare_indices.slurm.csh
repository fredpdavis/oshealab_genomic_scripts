#!/bin/csh -v
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-core=1
#SBATCH --mem=50g
#SBATCH -o slurm_out/prepare_indices.%A.out
#SBATCH -e slurm_out/prepare_indices.%A.err
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:200

# Load programs
module load kallisto/0.44.0
module load STAR/2.5.4a
module load ucsc

set BASEDIR="/data/davisfp/projects/cytokineX"
set NUMCPU=8

## SPECIFY INFORMATION ABOUT ASSEMBLy 
set ENSEMBL_GENOMEVER="GRCm38"
set ENSEMBL_RELEASE="94"
set EGR="${ENSEMBL_GENOMEVER}.${ENSEMBL_RELEASE}"
set GENOMEDIR="$BASEDIR/$EGR"
set GENOMEFA="$GENOMEDIR/genome.$EGR.fa"
set TRANSCRIPTFA="$GENOMEDIR/Mus_musculus.${ENSEMBL_GENOMEVER}.cdna.all.fa.gz"

set KALLISTODIR="$BASEDIR/kallisto_files.$EGR"
set KALLISTO_INDEX_NAME="$EGR.kalisto_index"

set STARDIR="$BASEDIR/star_files.$EGR"
set STAR_INDEX_NAME="$EGR.star_index"

set PICARDDIR="$BASEDIR/picard_files.$EGR"
set PICARD_REF_FLAT="refFlat.$EGR.txt"

set ORIG_GTF_FN="$GENOMEDIR/Mus_musculus.$EGR.gtf.gz"
set UNCOMPR_GTF_FN="Mus_musculus.$EGR.gtf"

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`

echo "# slurm run started on $curhost at $curtime"

## If not already here, download genome assembly and gene sets from ENSEMBL
if (! -s $GENOMEDIR ) then
   set curtime=`date`
   echo "Retrieving genome sequences and gene information ($curtime)"
   mkdir -p $GENOMEDIR
   cd $GENOMEDIR
   #1. GTF file
   wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus/Mus_musculus.$EGR.gtf.gz
   wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus/Mus_musculus.$EGR.chr_patch_hapl_scaff.gtf.gz
   
   #2. Transcript sequences
   wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/cdna/Mus_musculus.${ENSEMBL_GENOMEVER}.cdna.all.fa.gz
   
   #2b. noncoding RNA sequences
   wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/ncrna/Mus_musculus.${ENSEMBL_GENOMEVER}.ncrna.fa.gz

   #3. Genome sequences
   foreach chrom (`seq 1 19` 'X' 'Y' 'MT') 
      set CURCHROMSEQF="Mus_musculus.$ENSEMBL_GENOMEVER.dna_sm.chromosome.$CHROM.fa.gz"
      wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna/$CURCHROMSEQF
      zcat $CURCHROMSEQF >> $GENOMEFA
   end
   wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna/Mus_musculus.$ENSEMBL_GENOMEVER.dna_sm.nonchromosomal.fa.gz

   cd $curdir

endif

## If not already prepared, make KALLISTO index
if (! -s $KALLISTODIR) then
   set curtime=`date`
   echo "KALLISTO: preparing index"
   mkdir -p $KALLISTODIR
   cd $KALLISTODIR

   kallisto index -i $KALLISTO_INDEX $TRANSCRIPTFA
   zcat $ORIG_GTF_FN | grep '	transcript	' | sed 's/.*ENSMUSG/ENSMUSG/' | sed 's/".*ENSMUST/ ENSMUST/g' | sed 's/".*gene_name "/ /' | sed 's/".*//' | sort | uniq > tx_gene_name_map.txt

   cd $curdir
endif


## If not already prepared, make STAR index
if (! -s $STARDIR) then
   set curtime=`date`
   echo "STAR: preparing index"
   mkdir -p $STARDIR
   cd $STARDIR

   zcat $ORIG_GTF_FN > $UNCOMPR_GTF_FN
   STAR --runMode genomeGenerate --genomeDir $STAR_INDEX_NAME --genomeFastaFiles $GENOMEFA --runThreadN $NUMCPU --sjdbGTFfile $UNCOMPR_GTF_FN

   cd $curdir
endif


## If not already prepared, make gene structure file for PICARD
if (! -s $PICARDDIR) then
   set curtime=`date`
   echo "PICARD: preparing refflat"
   mkdir -p $PICARDDIR
   cd $PICARDDIR
   gtfToGenePred -genePredExt -geneNameAsName2 $UNCOMPR_GTF_FN refFlat.tmp.txt
   paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) | awk 'BEGIN{FS="\t"; OFS="\t"} {if ($1 == "") {$1 = $2} print $0}' > $PICARD_REF_FLAT
   gzip $PICARD_REF_FLAT

   cd $curdir
endif

set curtime=`date`
echo "# slurm run finished on $curhost at $curtime"
