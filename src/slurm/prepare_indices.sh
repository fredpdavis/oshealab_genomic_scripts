#!/bin/csh -v
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-core=1
#SBATCH --mem=50g
#SBATCH -o slurm_out/prepare_indices.%A.out
#SBATCH -e slurm_out/prepare_indices.%A.err
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:200

module load kallisto/0.44.0
module load STAR/2.5.4a

set BASEDIR=""

## If not already here, download genome assembly and gene sets from ENSEMBL
## Genome assembly: GRCm38.p4 -- ENSEMBL Release 82

set ENSEMBL_GENOMVER="GRCm38"
set ENSEMBL_RELEASE="82"
set GENOMEDIR="$BASEDIR/${ENSEMBL_GENOMEVER}.${ENSEMBL_RELEASE}"
set KALLISTODIR="$BASEDIR/kallisto_files"
set STARDIR="$BASEDIR/star_files"
set KALLISTO_INDEX_NAME="${ENSEMBL_GENOMEVER}.${ENSEMBL_RELEASE}.kalisto_index"
set STAR_INDEX_NAME="${ENSEMBL_GENOMEVER}.${ENSEMBL_RELEASE}.star_index"

if (! -s $GENOMEDIR ) then
   mkdir -p $BASEDIR/GRCm38.82
   #1. GTF file
   wget ftp://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/Mus_musculus.GRCm38.82.gtf.gz
   wget ftp://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/Mus_musculus.GRCm38.82.chr_patch_hapl_scaff.gtf.gz
   
   #2. Transcript sequences
   wget ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
   
   #2b. noncoding RNA sequences
   wget ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz
endif

if (! -s $KALLISTODIR) then
   kallisto index -i $KALLISTO_INDEX ../GRCm38.ENSEMBL82/Mus_musculus.GRCm38.cdna.all.fa.gz 

   cat /fdb/cellranger/refdata-cellranger-1.2.0/mm10/genes/genes.gtf | grep '	transcript	' | sed 's/.*ENSMUSG/ENSMUSG/' | sed 's/".*ENSMUST/ ENSMUST/g' | sed 's/".*gene_name "/ /' | sed 's/".*//' | sort | uniq > tx_gene_name_map.txt
endif


if (! -s $STARDIR) then
   STAR --runMode genomeGenerate --genomeDir $STAR_INDEX_NAME --genomeFastaFiles GRCm38.82.ERCC.fa --runThreadN 16 --sjdbGTFfile GRCm38.82.ERCC.gtf

endif


3. Genome sequences
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.1.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.2.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.3.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.4.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.5.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.6.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.7.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.8.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.9.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.10.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.11.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.12.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.13.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.14.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.15.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.16.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.17.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.18.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.19.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.X.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.Y.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.MT.fa.gz
ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.nonchromosomal.fa.gz
