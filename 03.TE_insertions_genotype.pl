#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;


#--Usage-----------------------------------------

my $usage=<<USAGE;

          ****** Genotype TE insertions using short reads ******

          perl $0   -Fasta  <ref.referenceTEinsertions_and_flanking1kb.fasta>  -leftRead <Sam1_1.fq.gz>  -rightRead <Sam1_2.fq.gz> -samId <Sam1>  -output <Sam1.refereceTEinsertion>  -script <the path to scripts>  -threads <threads>

          -Fasta	[required] the TE insertion and flanking sequences in fasta format
          -leftRead	[required] left read 
	  -rightRead	[required] right read
	  -samId	[required] sample name i.e. Sam1
	  -output	[required] TE genotype results
          -script	[required] the path to scripts
          -threads	[optional] threads  default: 6 cores

          Author: Xu Cai
          Bug report: caixu0518\@163.com          


USAGE

my ($in0, $in1, $in2, $in3, $out, $script, $threads);

GetOptions(
        "Fasta:s"        =>\$in0,
        "leftRead:s"     =>\$in1,
        "rightRead:s"    =>\$in2,
        "samId:s"        =>\$in3,
        "output:s"       =>\$out,
        "script:s"       =>\$script,
        "threads:s"      =>\$threads,

);

die $usage if (!defined $in0 || !defined $in1 || !defined $in2 || !defined $in3 || !defined $out || !defined $script);


my $flankingLen = 1000;
my $errorLen = 10;

   $threads = 6, if(not defined $threads);


my $leftread  = $in1;
my $rightread = $in2;

my $cmdString = "NA";
   if(not -e "$in0.amb"){
      $cmdString = "$script/bwa index $in0";
      print STDERR (localtime) . ": CMD: $cmdString\n";
      system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
   }

   $cmdString = "$script/bwa mem  -t  $threads  -T 20  -Y   $in0  $leftread  $rightread  >  $in3.sam";
   print STDERR (localtime) . ": CMD: $cmdString\n";
   system("$cmdString") == 0 or die "failed to execute: $cmdString\n";


my $samFile = $in3.".sam";


my %id2Len = ();
   &readFasta($in0, \%id2Len);

my $Pos3readId;
open IN0, $samFile;

while(<IN0>){
  next, if(/^@/);
  chomp;
  my @temp = split(/\t/, $_);
     
     my ($startPos, $endPos) = ($temp[3], $temp[3]);

     my $flag = 0;
     my @m = ();
    
     if($temp[5] =~ /(\d+)M/g || $temp[5] =~ /(\d+)D/g){
        $endPos = $startPos + $1;
        my $tempLen = $1;
 
        if($tempLen >= 20){
            my $tmp = $tempLen."M";
            if($temp[5] =~ /$tmp/){
               if(not exists $Pos3readId ->{$temp[0]}){
                  $Pos3readId ->{$temp[0]} ->{$temp[2]} = "$startPos:$endPos";
               }
               else{
                  $Pos3readId ->{$temp[0]} ->{$temp[2]} .= "\t"."$startPos:$endPos";
               }
            }
        }
        $startPos = $endPos+1;

     }
}
close IN0;


open OUT, ">$in3.readMapped.regions.Lst";
for my $key1(keys %{$Pos3readId}){
    for my $key2(keys %{$Pos3readId ->{$key1}}){
        my $seqLen = $id2Len{$key2};
        my ($break1, $break2) = ($flankingLen, $seqLen-$flankingLen);
          
        my @info = split("\t", $Pos3readId ->{$key1} ->{$key2});
        my @labels = ();

        my $flag = 0;
        for my $each(@info){
            next, if($each !~ /:/);
            my @pos = split(/:/, $each);
               @pos = sort {$a<=>$b } @pos;
            my $matchedLen = $pos[1]-$pos[0]+1;
            my $value = "NA";
            my $aLen;  
            my $bLen;
            my $cLen;  

            if($pos[0] <= $break1){

               if($pos[1] <= $break1){
                  $value = "aa";
               }
               if($pos[1] > $break1 && $pos[1] <= $break2){
                  $aLen = $break1-$pos[0]+1;                   
                  $bLen = $pos[1] - $break1;
                   
                  if($aLen >= $errorLen && $bLen >= $errorLen){
                     $value = "ab";
                  }elsif($aLen >= $errorLen && $bLen < $errorLen){
                     $value = "aa";
                  }
                  elsif($aLen <= $errorLen && $bLen >= $errorLen){
                     $value = "bb";
                  }
                  else{
                     $value = "NA";
                  }
                   
               }
               if($pos[1] > $break2){
                  $aLen = $break1-$pos[0]+1;
                  $bLen = $break2 - ($break1+1) + 1; 
                  $cLen = $pos[1] - ($break2+1) + 1;
                  
                  if($aLen >= $errorLen && $cLen >= $errorLen){
                     $value = "ac";
                  }elsif($aLen < $errorLen  && $cLen >= $errorLen){
                     $value = "bc";
                  }elsif($aLen >= $errorLen && $cLen < $errorLen){
                     $value = "ab";
                  }elsif($aLen < $errorLen  && $cLen < $errorLen){
                     $value = "bb"; 
                  }
                  else{
                     $value = "NA";
                  } 
               }
            }
            
	    if($pos[0] > $break1  && $pos[0] < $break2){
               if($pos[1] >  $break1 &&  $pos[1] <= $break2){
                  $value = "bb";
               }
               if($pos[1] > $break2){
                  $bLen = $break2 - $pos[0] + 1;
                  $cLen = $pos[1]-($break2+1) +1;
                  
                  if($bLen >= $errorLen && $cLen >= $errorLen){
                     $value = "bc";
                  } 
                  elsif($bLen >= $errorLen && $cLen < $errorLen){
                     $value = "bb";
                  }
                  elsif($bLen < $errorLen && $cLen >= $errorLen){
                     $value = "cc";
                  }
                  else{
                     $value = "NA";
                  }
               }
            } 

            if($pos[0] >= $break2){
               if($pos[1] > $break2){
                  $value = "cc";
               }
            }
            push(@labels, $value); 
        }
        my $string = join("", @labels);
        my $a = $string;
        my $b = $string;
        my $c = $string;
        my $flagA = 0;
           $a =~ s/a//g;
           $b =~ s/b//g;
           $c =~ s/c//g;
           $flagA = 1,   if($a eq "");
           $flagA = 1,   if($b eq "");
           $flagA = 1,   if($c eq "");
           print OUT join("\t", $key1, $key2, $seqLen, @info, @labels), "\n", if($flagA == 0);
    }
}
close OUT;
system("perl  $script/code_genotype.pl  $in3.readMapped.regions.Lst  $out");
system("rm -rf $in3.readMapped.regions.Lst");

#system("gzip $in3.readMapped.regions.Lst");

`rm  $samFile`, if(-e  $samFile);

##--------------------subs------------------------------
sub ArraySum {
    
    my ($array, $sum) = @_;

    for my $key1(@{$array}){
        $$sum += $key1;
    }

}

sub readFasta {

  my ($in,$id2Len) = @_;
  my $id2seq;
  open(my $SFR,$in);

  my $id;
  while($_=<$SFR>) {
    if(/^>([^\s^\n]+)\s*\n*/) {
      $id = $1;
      $id2seq->{$id} = "";
    }
    else {
      chomp;
      $id2seq->{$id} .= $_;
    }
  }
  close($SFR);

  for my $key(keys %{$id2seq}){
      $id2Len ->{$key} = length($id2seq ->{$key});
  }

}



