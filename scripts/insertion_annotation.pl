#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- Merged.SV.INS.gz.withMorethan80_TE_cov
my $in1 = $ARGV[1]; ##- Brapa_Chiifu_v3.1_genes.gff3
my $flankingLen = 2000;
my $out = $in0.".anno";


my $genefeature = $in1.".bed";
   &gff2bed($in1, $genefeature);

   &output($genefeature, $in0, $out);
 
   `rm $genefeature`;


sub output {

    my ($bedFile, $svbed, $output) = @_;

    my $gene;
    my $Upstream;
    my $Downstream;
    my $CDS;

    open IN2, $bedFile;
    while(<IN2>){
      chomp;
      my @temp = split(/\t/, $_);
      for(my $m=$temp[1]; $m<=$temp[2]; $m++){
          $gene ->{$temp[0]} ->{$m} = $temp[5],   if($temp[4] eq "gene");
          $Upstream ->{$temp[0]} ->{$m} = $temp[5],   if($temp[4] eq "Upstream"); 
          $Downstream ->{$temp[0]} ->{$m} = $temp[5],   if($temp[4] eq "Downstream");
          $CDS ->{$temp[0]} ->{$m} = $temp[5],   if($temp[4] eq "CDS");
      }
    }
    close IN2;

    my %results = ();
    open IN3, $svbed;
    while(<IN3>){
      chomp;
      my @temp = split(/\t/, $_);
      my $upPre = "-";     
      my $downPre = "-"; 
      my $genePre = "-";
      my $cdsPre = "-";
         $upPre = $Upstream ->{$temp[0]} ->{$temp[1]},   if(exists $Upstream ->{$temp[0]} ->{$temp[1]}); 
         $downPre = $Downstream ->{$temp[0]} ->{$temp[1]},  if(exists $Downstream ->{$temp[0]} ->{$temp[1]});        
         $genePre = $gene ->{$temp[0]} ->{$temp[1]},  if(exists $gene ->{$temp[0]} ->{$temp[1]});
         $cdsPre = $CDS ->{$temp[0]} ->{$temp[1]},  if(exists $CDS ->{$temp[0]} ->{$temp[1]});
         my $flag = 0;
            $flag += 1, if($upPre ne "-");
            $flag += 1, if($downPre ne "-");
            $flag += 1, if($genePre ne "-");
            $flag += 1, if($cdsPre ne "-"); 
            if($flag > 0){
               my @values = ($temp[0], $temp[1], $temp[2], $temp[3], $temp[4], $temp[5], $upPre, $downPre, $genePre, $cdsPre);
               $results{$temp[0]}{$temp[1]}{$temp[2]} = \@values;
            }
    }
    close OUT;
    close IN3;

    open OUT, ">$output";
    print OUT join("\t", "Chr", "Start", "End", "Type", "SVlen", "Upstream", "Downstream", "Gene", "CDS"), "\n";
    for my $key1(sort keys %results){
        for my $key2(sort {$a<=>$b} keys %{$results{$key1}}){
            for my $key3(sort {$a<=>$b} keys %{$results{$key1}{$key2}}){
                print OUT join("\t", @{$results{$key1}{$key2}{$key3}}), "\n";
            }
        }
    }
    close OUT;  

}


sub gff2bed {

    my ($gff3, $bedOut) = @_;

    my $geneid;
    open OUT, ">$bedOut";
    open IN0, $gff3;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
      if($temp[2] eq "gene"){
         $geneid = $temp[8];
         $geneid =~ s/ID=//;
         $geneid =~ s/;//;
         print OUT join("\t", $temp[0], $temp[3], $temp[4],  $temp[6], $temp[2], $geneid), "\n", if($temp[2] eq "gene");
         if($temp[6] eq "+"){
            print OUT join("\t", $temp[0], $temp[3]-$flankingLen, $temp[3]-1, $temp[6], "Upstream", $geneid), "\n";
            print OUT join("\t", $temp[0], $temp[4]+1, $temp[4]+$flankingLen, $temp[6], "Downstream", $geneid), "\n";         
         }
         if($temp[6] eq "-"){
            print OUT join("\t", $temp[0], $temp[4]+1, $temp[4]+$flankingLen, $temp[6], "Upstream", $geneid), "\n";
            print OUT join("\t", $temp[0], $temp[3]-$flankingLen, $temp[3]-1, $temp[6], "Downstream", $geneid), "\n"; 
         }
      }
      else{
         print OUT join("\t", $temp[0], $temp[3], $temp[4], $temp[6], $temp[2], $geneid), "\n", if($temp[2] eq "mRNA" || $temp[2] eq "CDS");
      }
    }
    close IN0;
    close OUT;

}



