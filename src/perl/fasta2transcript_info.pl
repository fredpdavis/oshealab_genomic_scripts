#!/usr/local/bin/perl
#fpd 150518_2144
# Purpose: extract transcript info from BDGP (+ERCC + INTACT) GTF File

use strict;
use warnings;
main() ;

sub main {

   my $usage = __FILE__." [biotype] < FASTA > map.txt" ;

   my $specs = {} ;
   if ($#ARGV >= 0) {$specs->{biotype} = $ARGV[0];}

   print join("\t", qw/transcript_id gene_id transcript_name gene_name gene_biotype chr start end/)."\n";
   while (my $line = <STDIN>) {
      if ($line !~ /^\>/) {next;}
      chomp $line;
      $line =~ s/\>//;
#      print "in on $line\n";

#>ENSMUST00000177564.1 cdna chromosome:GRCm38:14:54122226:54122241:1 gene:ENSMUSG00000096176.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:Trdd2 description:T cell receptor delta diversity 2 [Source:MGI Symbol;Acc:MGI:4439546]

      my ($txid, $geneid, $txname, $genename, $biotype, $chr, $start, $end) ;
      if ($line =~ /^ENS/) {
         $txid = $line ;
         $txid =~ s/[ \t].*// ;
         $txname = $txid ;

         ($geneid) = ($line =~ /gene:(ENS[^ ]*)/) ;
         ($genename) = ($line =~ /gene_symbol:([^ ]*)/) ;
         ($biotype) = ($line =~ /gene_biotype:([^ ]*)/) ;

         my @t = split(' ', $line) ;
         (undef, undef, $chr, $start, $end,undef) = split(':', $t[2]) ;
#         die "searching in $t[2] found $chr $start $end\n";


      } else {
         $txid = $line ;
         $txid =~ s/[ \t].*// ;
         ($txname, $geneid, $genename) = ($txid, $txid, $txid) ;
         $biotype = "unk" ;
         $chr = '';
         $start = '';
         $end = '' ;
      }

      if (exists $specs->{biotype}) {
         $biotype = $specs->{biotype} ;
      }

      print join("\t", $txid, $geneid, $txname, $genename, $biotype, $chr, $start, $end)."\n";
   }

}

#3R	FlyBase	transcript	722370	722621	.	-	.	gene_id "FBgn0085804"; gene_version "1"; transcript_id "FBtr0114258"; transcript_version "1"; gene_name "CR41571"; gene_source "FlyBase"; gene_biotype "pseudogene"; transcript_name "CR41571-RA"; transcript_source "FlyBase"; transcript_biotype "pseudogene";
#3R	FlyBase	transcript	835381	2503907	.	+	.	gene_id "FBgn0267431"; gene_version "1"; transcript_id "FBtr0346770"; transcript_version "1"; gene_name "CG45784"; gene_source "FlyBase"; gene_biotype "protein_coding"; transcript_name "CG45784-RA"; transcript_source "FlyBase"; transcript_biotype "protein_coding";

#>ENSMUST00000196221.1 cdna chromosome:GRCm38:14:54113468:54113476:1 gene:ENSMUSG00000096749.2 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:Trdd1 description:T cell receptor delta diversity 1 [Source:MGI Symbol;Acc:MGI:4439547]
#>ENSMUST00000179664.1 cdna chromosome:GRCm38:14:54113468:54113478:1 gene:ENSMUSG00000096749.2 gene_biotype:TR_D_gene transcript_biotype:processed_transcript gene_symbol:Trdd1 description:T cell receptor delta diversity 1 [Source:MGI Symbol;Acc:MGI:4439547]
#>ENSMUST00000178862.1 cdna chromosome:GRCm38:6:41542163:41542176:1 gene:ENSMUSG00000094569.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:Trbd2 description:T cell receptor beta, D region 2 [Source:MGI Symbol;Acc:MGI:4439727]
#>ENSMUST00000178537.1 cdna chromosome:GRCm38:6:41533201:41533212:1 gene:ENSMUSG00000095668.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:Trbd1 description:T cell receptor beta, D region 1 [Source:MGI Symbol;Acc:MGI:4439571]
#>ENSMUST00000179520.1 cdna chromosome:GRCm38:12:113430528:113430538:-1 gene:ENSMUSG00000094028.1 gene_biotype:IG_D_gene transcript_biotype:IG_D_gene gene_symbol:Ighd4-1 description:immunoglobulin heavy diversity 4-1 [Source:MGI Symbol;Acc:MGI:4439801]
#>ENSMUST00000179883.1 cdna chromosome:GRCm38:12:113448214:113448229:-1 gene:ENSMUSG00000094552.1 gene_biotype:IG_D_gene transcript_biotype:IG_D_gene gene_symbol:Ighd3-2 description:immunoglobulin heavy diversity 3-2 [Source:MGI Symbol;Acc:MGI:4439707]
#>ENSMUST00000195858.1 cdna chromosome:GRCm38:12:113449588:113449597:-1 gene:ENSMUSG00000096420.2 gene_biotype:IG_D_gene transcript_biotype:IG_D_gene gene_symbol:Ighd5-6 description:immunoglobulin heavy diversity 5-6 [Source:MGI Symbol;Acc:MGI:4937234]
#>ENSMUST00000179932.1 cdna chromosome:GRCm38:12:113449588:113449599:-1 gene:ENSMUSG00000096420.2 gene_biotype:IG_D_gene transcript_biotype:processed_transcript gene_symbol:Ighd5-6 description:immunoglobulin heavy diversity 5-6 [Source:MGI Symbol;Acc:MGI:4937234]
#>ENSMUST00000180001.1 cdna chromosome:GRCm38:12:113450851:113450867:-1 gene:ENSMUSG00000095656.1 gene_biotype:IG_D_gene transcript_biotype:IG_D_gene gene_symbol:Ighd2-8 description:immunoglobulin heavy diversity 2-8 [Source:MGI Symbol;Acc:MGI:4439706]
