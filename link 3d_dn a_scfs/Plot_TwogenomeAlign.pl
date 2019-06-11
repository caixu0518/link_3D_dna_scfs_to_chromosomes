#!/usr/bin/perl -w
use strict;
use lib "/home/caix/perl5/lib/perl5";
use SVG;
require("/home/caix/bin/all-sub.pl");

my $in0 = $ARGV[0]; # Ref Chr
my $in1 = $ARGV[1]; # Query Chr
my $index = $ARGV[2]; # for example: JZS-0212

my $pathNucmer = "/data/mg1/caix/src/biosoft/mummer-4.0.0beta2";     # it depends__
my $pathChr2SVG = "/data/mg1/caix/scripts/draw/nucmer";    # it depends__

###---###---###---###
my $out0 = $in0.".Len";
   &ChrLen($in0, $out0);

my $out1 = $in1.".Len";
   &ChrLen($in1, $out1);

system("$pathNucmer/nucmer -c 4000 -g 500 -l 2000 -t 30 -p $index  $in0  $in1");

system("$pathNucmer/show-coords -r $index.delta > $index.coords");

my $sortedLoci = "$index.coords.loci.sort";
   &Coords_Loci("$index.coords", $sortedLoci);

system("perl $pathChr2SVG/Chr2ChrSVG.ppt.pl  $out0 $out1 $sortedLoci");  # X axis chromosome length file   y axis chromosome length file coords sorted file

###subs###
sub Coords_Loci {
  my ($coordF, $out) = @_;  #get loci and sort
  my ($temp) = $coordF.".Temp";
 
  open IN0, $coordF;
  open OUT0, ">$temp";
  my $flag=0;
  while (<IN0>){
    chomp;
    $flag++;
    if ($flag>=6){
      my @info=split;
      #my $ref=($info[0]+$info[1])/2;
      #my $query=($info[3]+$info[4])/2;
      print OUT0 "$info[11]\t$info[0]\t$info[1]\t$info[12]\t$info[3]\t$info[4]\n"; ##- ref s e query s e
    }
  }
  close IN0;
  close OUT0;
  system("sort -k1,1 -k2,2n $temp > $out");  
  system("rm -rf $temp");

}
