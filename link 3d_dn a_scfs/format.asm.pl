#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- fasta
my $in1 = $ARGV[1]; ##- matched list (e.g. 1	A01) 
my $out = "final.asm.chr.fasta";


##------------------------------------------------------------------------------
##- read matched file
my %match = ();
   &read_match($in1, \%match);

##- read Seq
my %id2seq = ();
my @order = ();
   &readFasta($in0, \%id2seq, \%match, \@order);
   #print ">", $_, "\n", $id2seq{$_}, "\n", for (keys %id2seq);

##- generate Output
my $tmpout = $out.".tmp";
    &output(\%id2seq, \@order, $tmpout);
  
##- format fasta assembly
    `fold $tmpout > $out`;

##- clean
    `rm $tmpout`;


##----------------------------------------------------------------------
sub output {

    my ($id2seq, $order, $out) = @_;
    
    my $allnum = scalar(@{$order}) - 1;
    print "Number of total sequences:", $allnum + 1, "\n";
    my $scfCount = 0;
    open OUT0, ">$out";
    for(my $m=0; $m<=$allnum; $m+=1) {
        my $id = $order ->[$m];

        if($m<=9){  ##- it depends_  9 means chromosomes number -1 (eg. 10 for 9)
           print OUT0 ">", $id, "\n", $id2seq ->{$id}, "\n";
        }
        else{
           $scfCount += 1;
           my $numfmt;
           &formatNum($scfCount, \$numfmt);
           $numfmt = "Scaffold".$numfmt;
           print OUT0 ">", $numfmt, "\n", $id2seq ->{$id}, "\n";
        }
    }
    close OUT0;

}

sub read_match {
  
    my ($matchedFile, $match) = @_;
 
    open IN1, $matchedFile;
    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_);
      $match ->{$temp[0]} = $temp[1];
    }
    close IN1;

}

sub readFasta {

    my ($fasta, $id2seq, $match, $order) = @_;

    my (@Chr, @Scaffolds);
    my $id = 'NA';
    open IN0, $fasta;
    while(<IN0>){
      chomp;
      if(/^>(\S+)/){
         $id = $1;
         if(exists $match ->{$1}){
            $id = $match ->{$1};
            push(@Chr, $id);
         }
         else{
            push(@Scaffolds, $id);
         }
         $id2seq ->{$id} = "";
      }
      else{
         $id2seq ->{$id} .= $_;
      }
    }
    close IN0;

    ##- sort scaffolds
    my $scafffoldIndex;
     
    for my $eachScf(@Scaffolds){
        my $len = length($id2seq ->{$eachScf});
        $scafffoldIndex ->{$len} ->{$eachScf} = "Y";
    } 
  
    my @newScfOrder = ();
    for my $key1(sort {$b<=>$a} keys %{$scafffoldIndex}){
        for my $key2(sort keys %{$scafffoldIndex ->{$key1}}){
            push(@newScfOrder,  $key2); 
        }
    }
    push(@{$order}, sort(@Chr), @newScfOrder);
   
}

sub formatNum {
  my ($n, $num) = @_;

  $$num = "000".$n  if($n >= 1 && $n < 10);
  $$num = "00".$n  if($n >= 10 && $n < 100);
  $$num = "0".$n  if($n >= 100 && $n < 1000);
  $$num = $n  if($n >= 1000);

}

