#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- All sampleName
my $out = "TIP.array";

my @Sps = ();
open IN0, $in0;
while(<IN0>){
  chomp;
  my $id = $_;
  push(@Sps, $id);
}
close IN0;

my %index = ();
for my $key(@Sps){
  
    my $file = "$key.Non-refereceTEinsertion";  ##- it depends_ 

    if(-e $file){
       open IN1, $file;
       <IN1>;
       while(<IN1>){
         chomp;
         my @temp = split(/\t/, $_);
         $index{$temp[0]}{$key} = $temp[6];
       }
       close IN1;
    }
    else{
       die "Error: cannot find $file\n";
    } 
}

open OUT, ">$out";
print OUT join("\t", "TEindex", @Sps), "\n";
for my $key(sort keys %index){
    my @lineInfo = ();
    push(@lineInfo, $key);
    for my $each(@Sps){
        my $value = "NA";
           $value = $index{$key}{$each}, if(exists $index{$key}{$each}); 
           push(@lineInfo, $value);
    }
    print OUT join("\t", @lineInfo), "\n"; 
}
close OUT;
