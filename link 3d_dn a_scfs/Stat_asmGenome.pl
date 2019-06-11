#!/usr/bin/perl -w
# Xu Cai
use strict;

my $in0 = $ARGV[0]; ##- final assembly fasta
my $chrNum = $ARGV[1]; ## Number of chromosomes

`samtools faidx $in0`;

`cut -f 1-2 $in0.fai | sort -k2,2nr > $in0.sizes`;

my ($totalSize, $chrSize, $scfSize, $ratio);
my $num = 0;
open IN0, "$in0.sizes";
while(<IN0>){
  chomp;
  my @temp = split(/\t/, $_);
     $num += 1;
  if($num <= $chrNum){
     $chrSize += $temp[1];
  }
  else{
     $scfSize += $temp[1];
  }
  $totalSize += $temp[1];
}
close IN0;

print "Total bases of final assembly: $totalSize (bp)\n";
print "Total bases of anchored sequeces (chromosomes): $chrSize (bp)\n";
print "Total bases of un-anchored sequeces (scaffolds): $scfSize (bp)\n";
print "Anchored ratio of the final assembly: ", $chrSize/$totalSize, "\n";
