#!/usr/bin/perl

use warnings;
use strict;


my $in0 = $ARGV[0]; ##- all.querys.lst  one column, describe the query name
my $prefix = $ARGV[1];  ##- the ref name
my $INSout = $prefix.".merged.INS";
my $DELout = $prefix.".merged.DEL";


my %INS = ();
my %DEL = ();
my %index = ();

open IN0, $in0;
while(<IN0>){
  chomp;
  my @temp = split(/\t/, $_);
  my $spid = $temp[0];
  my $fileName = "$prefix-".$spid.".svs.bed.gz";

  open IN1, "gzip -dc $fileName |";
  <IN1>;
  while(<IN1>){
    chomp;
    my @temp = split(/\t/, $_);
    my @result = ($temp[0],$temp[1],$temp[2],$temp[3],$temp[4], $temp[13]);
    $index{$temp[0]}{$temp[1]}{$temp[2]}  = \@result;  

    if($temp[3] eq "insertion"){
       if(not exists $INS{$temp[0]}{$temp[1]}{$temp[2]}){
          $INS{$temp[0]}{$temp[1]}{$temp[2]} = $spid;
       }
       else{
          $INS{$temp[0]}{$temp[1]}{$temp[2]} .= ";".$spid;
       }
    }   
     
    if($temp[3] eq "deletion") {
       if(not exists $DEL{$temp[0]}{$temp[1]}{$temp[2]}){
          $DEL{$temp[0]}{$temp[1]}{$temp[2]} = $spid;
       }
       else{
          $DEL{$temp[0]}{$temp[1]}{$temp[2]} .= ";".$spid;
       }
    }
  }
  close IN1;
}
close IN0;

open OUTINS, ">$INSout";
for my $key1(sort keys %INS){
    for my $key2(sort {$a<=>$b} keys %{$INS{$key1}}){
        for my $key3(sort {$a<=>$b} keys  %{$INS{$key1}{$key2}}){
            my $value = $INS{$key1}{$key2}{$key3};
            my @info = @{$index{$key1}{$key2}{$key3}};
            print OUTINS join("\t", $info[0], $info[1],$info[2],$info[3],$info[4],$value, $info[5]), "\n";
        }
    }
}
close OUTINS;

open OUTDEL, ">$DELout";
for my $key1(sort keys %DEL){
    for my $key2(sort {$a<=>$b} keys %{$DEL{$key1}}){
        for my $key3(sort {$a<=>$b} keys  %{$DEL{$key1}{$key2}}){
            my $value = $DEL{$key1}{$key2}{$key3};
            my @info = @{$index{$key1}{$key2}{$key3}};
            print OUTDEL join("\t", $info[0], $info[1],$info[2],$info[3],$info[4],$value, $info[5]), "\n";
        }
    }
}
close OUTDEL;

system("gzip $INSout");
system("gzip $DELout");
