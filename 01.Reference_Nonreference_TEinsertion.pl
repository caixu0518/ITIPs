#!/sur/bin/perl

use warnings;
use strict;
use Getopt::Long;


#--Usage-----------------------------------------

my $usage=<<USAGE;

        ****** Identification of reference and non-reference TE insertions between different genomes ******

	Usage: perl $0  -query <query.info.lst>  -ref <reference.info.lst>  -TElib <EDTA.TElib.fa>  -bin <the path to smartie-sv>  -script <the path to scripts> 


	-query	[required] the query id and query genome files. Two columns (queryName queryGenomeFile). 
        -ref    [required] the reference information. Three columns (referenceName ReferenceGenomeFile ReferenceGff3)
        -TElib  [required] the species TE library
        -bin 	[required] the path to smartie-sv. i.e. /10t/caix/src/smartie-sv/bin
        -script [required] the path to perl scripts 

        Author: Xu Cai
        Bug report: caixu0518\@163.com


USAGE


my ($in0, $in1, $in2, $smartieSv, $script);

GetOptions(
        "query:s"   =>\$in0,
        "ref:s"     =>\$in1,
        "TElib:s"   =>\$in2,
        "bin:s"     =>\$smartieSv,
        "script:s"  =>\$script,
);
die $usage if (!defined $in0 || !defined $in1);






my @queryIds = ();
my @queryGenomes = ();
my $refId;
my $refGenome;
my $refGff3;
   &readQueryRef($in0, $in1, \@queryIds, \@queryGenomes, \$refId, \$refGenome, \$refGff3);

   &main();
 

sub main {

    my $cmdString = "NA";

    ##- step 1: Using smartie-sv to detect insertions and deletions
    for(my $m=0; $m<=$#queryIds; $m++){
        my $queryName   = $queryIds[$m];
        my $queryGenomeFile  = $queryGenomes[$m];
        $cmdString = "perl  $script/Run_smartie-sv.pl  $queryGenomeFile  $queryName  $refGenome   $refId  $smartieSv";
        print STDERR (localtime) . ": COMMAND: $cmdString\n";
        system("$cmdString") == 0 || die "failed: $cmdString\n";
    }
     
    ##- step 2: Merge insertions and deletions 
    $cmdString = "perl $script/Merge_insertions_deletions.pl  $in0  $refId"; 
    print STDERR (localtime) . ": COMMAND: $cmdString\n";
    system("$cmdString") == 0 || die "failed: $cmdString\n";

    ##- step 3: blast the insertions and deletions to the species TE library
    $cmdString = "perl $script/blast_INS_DEL_seq_to_TE_lib.pl  $refId.merged.DEL.gz  $in2";
    print STDERR (localtime) . ": COMMAND: $cmdString\n";
    system("$cmdString") == 0 || die "failed: $cmdString\n";
    
    $cmdString = "perl $script/blast_INS_DEL_seq_to_TE_lib.pl  $refId.merged.INS.gz  $in2";
    print STDERR (localtime) . ": COMMAND: $cmdString\n";
    system("$cmdString") == 0 || die "failed: $cmdString\n";

    ##- step 4: extract TE insertions in the genic regions. 2kb upstream and downstream of the gene body
    ##- annotate reference TE insertions
    $cmdString = "perl  $script/deletion_annotation.pl  $refId.merged.DEL.gz.withMorethan80_TE_cov  $refGff3";
    print STDERR (localtime) . ": COMMAND: $cmdString\n";
    system("$cmdString") == 0 || die "failed: $cmdString\n";

    ##- annotate non-reference TE insertions
    $cmdString = "perl  $script/insertion_annotation.pl  $refId.merged.INS.gz.withMorethan80_TE_cov  $refGff3";
    print STDERR (localtime) . ": COMMAND: $cmdString\n";
    system("$cmdString") == 0 || die "failed: $cmdString\n";
 
    system("mv $refId.merged.DEL.gz.withMorethan80_TE_cov.anno  reference_TE.insertions.xls");
    system("mv $refId.merged.INS.gz.withMorethan80_TE_cov.anno  non-reference_TE.insertions.xls");
    
    ##- clean
    system("rm $in2.*  *.svs.bed.gz");

}



sub readQueryRef {

    my ($query, $ref, $queryIds, $queryGenomes, $refId, $refGenome, $refGff3) = @_;

    open IN0Q, $query;
    while(<IN0Q>){
      chomp;
      my @temp = split(/\t/, $_);
      push(@{$queryIds}, $temp[0]);
      push(@{$queryGenomes}, $temp[1]);
    }
    close IN0Q;

    open IN1R, $ref;
    while(<IN1R>){
      chomp;
      my @temp = split(/\t/, $_);
      $$refId = $temp[0];
      $$refGenome = $temp[1];
      $$refGff3 = $temp[2];
    }
    close IN1R;

}

