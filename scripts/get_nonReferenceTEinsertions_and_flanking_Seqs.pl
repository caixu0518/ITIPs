#!/usr/bin/perl

use warnings;
use strict;

my $in0 = $ARGV[0]; ##- genome fasta
my $in1 = $ARGV[1]; ##- Merged.SV.INS.gz.withMorethan80_TE_cov.anno
my $in2 = $ARGV[2]; ##- Merged.SV.INS.gz.withMorethan80_TE_cov
my $out = $ARGV[3]; ##- refId  i.e. ref 

my $flankLen = 1000;
my $outFasta = $out.".Non-referenceTEinsertions_and_flanking1kb.fasta";
my $outInfo  = $out.".Non-referenceTEinsertions_and_flanking1kb.info";


my %chr2seq = ();
my %Chr2Len = ();
   &readFasta($in0, \%chr2seq, \%Chr2Len);

my %TEPosIndex = ();
   &readPos($in1, $in2, \%TEPosIndex);

   &output(\%TEPosIndex, \%chr2seq, \%Chr2Len, $outFasta, $outInfo);

sub output {

    my ($TEPosIndex, $chr2seq, $Chr2Len, $outFasta, $outInfo) = @_;

    open OUTFasta, ">$outFasta";
    open OUTinfo, ">$outInfo";
    my $count = 0;
    for my $chr(keys %{$TEPosIndex}){
        for my $start(keys %{$TEPosIndex ->{$chr}}){
            for my $end(keys %{$TEPosIndex ->{$chr} ->{$start}}){
                        $count += 1;
                        my @temp = @{$TEPosIndex ->{$chr} ->{$start} ->{$end}};
                        my $chrLen = $Chr2Len ->{$chr}; 
                        my $upstart = $start - $flankLen+1;
                           $upstart = 0, if($upstart < 0);
                        my $downStop = $end + $flankLen;
                           $downStop = $chrLen, if($downStop > $chrLen);
                        my $upSeq = substr($chr2seq ->{$chr}, $upstart, ($start-$upstart+1));
                        my $downSeq = substr($chr2seq ->{$chr}, $end+1, $downStop-($end+1)+1);
                        my $midSeq = pop(@temp);                       
 
                        my $typeShort = $temp[3];
                           $typeShort = "D", if($temp[3] eq "deletion");
                           $typeShort = "I", if($temp[3] eq "insertion");
                       
                        my $indexName = $typeShort.$out.$count;
                        my $seq = $upSeq.$midSeq.$downSeq;
                        print OUTFasta  ">$indexName\n";
                        print OUTFasta  $seq, "\n"; 
                        print OUTinfo   join("\t", $indexName, @temp), "\n";
            }
        }
    }
    close OUTFasta;
    close OUTinfo; 

}

sub readPos {

    my ($in0, $svFiles, $index) = @_;

    open IN1, $in0;
    <IN1>;
    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_);
      $index ->{$temp[0]} ->{$temp[1]} ->{$temp[2]} = \@temp;
    }
    close IN1;

    open IN2, $svFiles;
    while(<IN2>){
      chomp;
      my @temp = split(/\t/, $_);
      if(exists $index ->{$temp[0]} ->{$temp[1]} ->{$temp[2]}){
         my @value = @{$index ->{$temp[0]} ->{$temp[1]} ->{$temp[2]}};
         push(@value, $temp[-1]);
         $index ->{$temp[0]} ->{$temp[1]} ->{$temp[2]} = \@value;
      }
    }
    close IN2;

}


sub readFasta {

  my ($in,$id2seq, $che2Len) = @_;
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
      $che2Len ->{$key} = length($id2seq ->{$key});
  }

}






