#!/usr/bin/perl -w
# Xu Cai
use strict;
require("/home/caix/bin/all-sub.pl");

my $in0 = $ARGV[0]; ##- corrected fasta
my $in1 = $ARGV[1]; ##- *review.assembly
my $gapSize = 100;

   &main();
sub main {

  my %id2seq = ();
     &readFasta($in0, \%id2seq);
  ##---------------------------------------------
  my $gap;
  for(my $n =1; $n <= $gapSize; $n++){
      $gap .= "N";
  }

  ##----------------------------------------------
  my (%pos2Scf, %length, %genome);
  my $out = "$in0.linkToChr.fa";
     &readReviewAsm(\%id2seq, \%pos2Scf, \%length, \%genome, $gap, $in1, $out);
  

}

sub readReviewAsm {
   my ($id2seq, $pos2Scf, $length, $genome, $gap, $ReviewAsmFile, $out)  = @_;
   
   ##-----------index--------------------------------------------------------------
   open IN0, $ReviewAsmFile;
   my $num = 0;
   my $scfNum = 10;
   while(<IN0>){
      if(/^>/){ ##- process pos to scaffold
         chomp;
         my @temp = split(/\s+/, $_);
         my $id = $temp[0];
            $id =~ s/^>//;
            $pos2Scf ->{$temp[1]} = $id, if(not exists $pos2Scf ->{$temp[1]});
            $length ->{$temp[1]} = $temp[3], if(not exists $length ->{$temp[1]});
      }
      else{
            $num += 1;
            chomp;
            if($num <= 10){ ##- chrmosome numbers
               $genome ->{$num} = $_, if(not exists $genome ->{$num});  
            }
            else{
               my @scfs = split(/\s+/, $_);
               for my $eachScf(@scfs){
                   $scfNum += 1;
                   $genome ->{$scfNum} = $eachScf;
               }
            }
      } 
   } 
   close IN0;

   ##--------------------------------------------------------------------------------
   open OUT0, ">$out";
   for my $key(sort {$a<=>$b} keys %{$genome}){
     my @element = split(/\s+/, $genome ->{$key});
     if(scalar(@element) >= 2){
        my @seq = ();
        my $seq;
        for my $pos(@element){
               $seq = $id2seq ->{$pos2Scf ->{abs($pos)}};
               if($pos < 0){
                  $seq =~ tr/ATCGatcg/TAGCtagc/;
                  $seq = reverse($seq);
               }
               push(@seq, $seq);
        }
        my $linkedSeq = join($gap, @seq);
        print OUT0 ">", $key, "\n", $linkedSeq, "\n";
     }
     else{
        my $pos = $element[0];
        my $seq = $id2seq ->{$pos2Scf ->{abs($pos)}};
        print OUT0 ">", $key, "\n", $seq, "\n";
     }
   }
   close OUT0;
   
}

