#!/use.binperl -w
use strict;
##- Xu Cai

my $in0 = $ARGV[0]; ##- primers   include Up and down sequences 
my $ratio = 0.5;
my $cutOffLen = 15;

   &process();
sub process {

   my $UpDownSeq = "UpDown.fa";
      &getFasta($in0, $cutOffLen, $UpDownSeq);

   my %id2seq = ();
      &indexFasta($UpDownSeq, \%id2seq);

   my %mark = ();
      &caculateRatio(\%id2seq, $ratio, \%mark);

   ##- remove low complexity sequences
   my $out =  "$in0.rm.xls";
      &rm($in0, \%mark, $out); 

}

sub rm {

    my ($primers, $mark, $out) = @_;

    my %index = ();
    for(keys %{$mark}){
       my @temp = split(/_/, $_);
       $index{$temp[0]}{$temp[1]} = "Y";
    }

    open IN1, $primers;
    open OUT0, ">$out";
    while(<IN1>){
      chomp;
      my $flag = "PASS";
      my @temp = split(/\s+/, $_);
         $flag = "Fail", if(exists $index{$temp[0]}{$temp[1]});
         print OUT0 join("\t", @temp, $flag), "\n";
    }
    close IN1;
    close OUT0;
}


sub getFasta {
     
    my ($primers, $cutOffLen, $outFasta) = @_;

    open IN0, $primers;
    open OUT0, ">$outFasta";
    while(<IN0>){
      chomp;
      my @temp = split(/\s+/, $_);
      my ($id5, $id3) = ("$temp[0]_$temp[1]_5", "$temp[0]_$temp[1]_3");
         $temp[-2] = reverse($temp[-2]);
      my $seq5 = substr($temp[-2], 0, $cutOffLen);
      my $seq3 = substr($temp[-1], 0, $cutOffLen);
      print OUT0 ">$id5\n$seq5\n";
      print OUT0 ">$id3\n$seq3\n";
    }
    close IN0;
    close OUT0;

}


sub caculateRatio {
   
    my ($id2seq, $ratio, $mark) = @_;
        
    
    for my $key1(sort keys %{$id2seq}){
        
        ##- one base repeat
        my $value1 = $id2seq ->{$key1}; $value1 =~ s/(\w{1})/$1\t/g;
        my @value1 = split(/\t/, $value1);
        my %count1 = ();
           $count1{$_} += 1, for(@value1);
           for(keys %count1){
               $mark ->{$key1} = "Y"  if($count1{$_}/scalar(@value1) >= $ratio);
           }

        ##- two bases repeat
        my $value2 = $id2seq ->{$key1}; $value2 =~ s/(\w{2})/$1\t/g;
        my @value2 = split(/\t/, $value2);
        my %count2 = ();
           $count2{$_} += 1, for(@value2);
           for(keys %count2){
               $mark ->{$key1} = "Y"  if($count2{$_}/scalar(@value2) >= $ratio);
           }
            
        ##- three bases repeat 
        my $value3 = $id2seq ->{$key1}; $value3 =~ s/(\w{3})/$1\t/g;
        my @value3 = split(/\t/, $value3);        
        my %count3 = ();
           $count3{$_} += 1, for(@value3);
           for(keys %count3){
               $mark ->{$key1} = "Y"  if($count3{$_}/scalar(@value3) >= $ratio);
           }

        ##- four bases repeat
        my $value4 = $id2seq ->{$key1}; $value4 =~ s/(\w{4})/$1\t/g;
        my @value4 = split(/\t/, $value4);
        my %count4 = ();
           $count4{$_} += 1, for(@value4);
           for(keys %count4){
               $mark ->{$key1} = "Y"  if($count4{$_}/scalar(@value4) >= $ratio);
           }
        
        ##- five bases repeat
        my $value5 = $id2seq ->{$key1}; $value5 =~ s/(\w{5})/$1\t/g;
        my @value5 = split(/\t/, $value5);
        my %count5 = ();
           $count5{$_} += 1, for(@value5);
           for(keys %count5){
               $mark ->{$key1} = "Y"  if($count5{$_}/scalar(@value5) >= $ratio);
           }  
    }
    
}


sub indexFasta {

    my ($fasta, $index) = @_;

    open IN0, $fasta;
    while(<IN0>){
       chomp;
       if(/^>(\S+)/){
          $index ->{$1} = "";
       }
       else{
          $index ->{$1} .= $_;
       }
    }
    close IN0;

}


sub Run_CMD {
   
    my ($cmdString) = @_;
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n"; 

}


