#!/usr/bin/perl

use warnings;
use strict;

my $in0 = $ARGV[0]; ##- genome fasta
my $in1 = $ARGV[1]; ##- reference_TE.insertions.xls
my $out = $ARGV[2]; ##- reference id

   $out = "Final",  if(not defined $out);

my $flankLen = 1000;
my $outFasta = $out.".referenceTEinsertions_and_flanking1kb.fasta";
my $outInfo  = $out.".referenceTEinsertions_and_flanking1kb.info";

my %chr2seq = ();
my %Chr2Len = ();
   &readFasta($in0, \%chr2seq, \%Chr2Len);

my %TEPosIndex = ();
   &readPos($in1, \%TEPosIndex);


   &output(\%chr2seq, \%Chr2Len, \%TEPosIndex, $outFasta, $outInfo);

sub output {

    my ($chr2seq, $Chr2Len, $TEPosIndex, $outFasta, $outInfo) = @_;

    my $count = 0;
    open OUTFasta, ">$outFasta";
    open OUTinfo, ">$outInfo";
    for my $chr(keys %{$TEPosIndex}){
        for my $start(keys %{$TEPosIndex ->{$chr}}){
            for my $end(keys %{$TEPosIndex ->{$chr} ->{$start}}){
                my @temp = @{$TEPosIndex ->{$chr} ->{$start} ->{$end}};
                $count += 1;
                my $chrLen = $Chr2Len ->{$chr};
                my $upstart = $start - $flankLen;
                   $upstart = 0, if($upstart < 0);
                my $downStop = $end + $flankLen;
                   $downStop = $chrLen, if($downStop > $chrLen);

                my $typeShort = $temp[3];
                $typeShort = "D", if($temp[3] eq "deletion");
                $typeShort = "I", if($temp[3] eq "insertion"); 
                
                my $upSeq = substr($chr2seq ->{$chr}, $upstart, ($start-1-$upstart+1));
                my $midSeq = substr($chr2seq ->{$chr}, $start, ($end-$start+1));
                my $downSeq = substr($chr2seq ->{$chr}, $end+1, ($downStop-($end+1)+1));

                my $indexName = $typeShort.$out.$count;
                print OUTFasta  ">$indexName\n";
                print OUTFasta  $upSeq.$midSeq.$downSeq, "\n";
                print OUTinfo   join("\t", $indexName, @temp), "\n";
            }
        }
    }
    close OUTFasta;
    close OUTinfo;


}


sub readPos {

    my ($in0, $index) = @_;

    open IN1, $in0;
    <IN1>;
    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_);
      $index ->{$temp[0]} ->{$temp[1]} ->{$temp[2]} = \@temp;
    }
    close IN1;

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
