#!/usr/bin/perl -w
use strict;
use Cwd;

##-source activate py3.5

my $in0    = $ARGV[0]; ##- query fasta
my $in1    = $ARGV[1]; ##- query prefix
my $in2    = $ARGV[2]; ##- reference fasta
my $in3    = $ARGV[3]; ##- reference prefix
my $bin    = $ARGV[4]; ##- the path where the smartie-sv program is located   i.e.  /10t/caix/src/smartie-sv/bin  

my $pwd = getcwd;
my $prefix = $in1."_".$in3;
my $smartie_sv = $bin;
 
   ##- check if the program exists
   my $program = "$smartie_sv/sawriter";
   die "Error: cannot find: $program\n", if(not -e $program);
   $program = "$smartie_sv/blasr";  
   die "Error: cannot find: $program\n", if(not -e $program);
   $program = "$smartie_sv/printgaps";
   die "Error: cannot find: $program\n", if(not -e $program); 

   system("rm -rf $prefix"), if(-d $prefix);
   system("mkdir $prefix");
   chdir $prefix;
   system("ln -s ../$in0 .");
   system("ln -s ../$in2 .");

   `cp  $smartie_sv/config.sh  .`;
   `cp  $smartie_sv/Snakefile  .`;
   `cp $smartie_sv/config.json .`;

    my $contigs = $in0.".contigs.fa";
       &scf2contig($in0, $contigs);
    

##- modify config.json
    `sed -i -r 's/hg38/$in3/' config.json`;
    `sed -i -r 's/hg38-chr19_39.4_39.8.fa/$in2/'  config.json`;
    `sed -i -r 's/gorilla/$in1/'  config.json`;
    `sed -i -r 's/000943F_quiver_patched.fa_0.5_0.9.fa/$contigs/'  config.json`;

`$smartie_sv/sawriter $in2`, if(not -e "$in2.sa");
`snakemake -s Snakefile -w 50  -p -k -j 20`;

##- clean
`rm $in2.fai  $in2.sa`;
#`rm  -rf  mappings  unmappings `;

chdir $pwd;

my $svfileName = $in3."-".$in1.".svs.bed";
`cp    $pwd/$prefix/variants/$svfileName   .`;
`gzip  $svfileName`;


sub scf2contig {

    my ($scaf, $out) = @_;

    my $id2seq;
       open IN0, $scaf;
       while(<IN0>){
         chomp;
         if(/^>(\S+)/){
            $id2seq ->{$1} = "";
         }
         else{
            $id2seq ->{$1} .= $_;
         }
       }
       close IN0;

    open OUT0, ">$out";
    for my $key(sort keys %{$id2seq}){
        my @temp = split(/N+/, $id2seq ->{$key});
        for(my $m=0; $m<= $#temp; $m++){
            print OUT0 ">", $key."_".$m, "\n", $temp[$m], "\n";
        }
    }
    close OUT0;

}
