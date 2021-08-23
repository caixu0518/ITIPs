# ITIPs
## Identification TE insertions polymorphisms based on a pan-genome and large-scale resequencing data


# Introduction
Identification TE insertions polymorphisms based on a pan-genome and large-scale resequencing data
## Step 1：identification of reference insertions and deletions.

Using each of the 20 *B. rapa* genomes as the reference and aligned each of the left 19 genome sequneces to the reference genome and call reference insertions and deletions using the [smartie-sv](https://github.com/zeeev/smartie-sv) piepeline.

![image](https://github.com/caixu0518/ITIPs/blob/main/png/Step%201.png)

## Step 2: Determination of TE insertions polymorphisms (TIPs).

![image](https://github.com/caixu0518/ITIPs/blob/main/png/Step%202.png)

## Step 3: Genotype TE insertions in a large-scale *B. rapa* population using short reads.
![image](https://github.com/caixu0518/ITIPs/blob/main/png/Step%203.png)

# Installation
The pipeline ITIPs is installation-free but requires dependencies: [smartie-sv](https://github.com/zeeev/smartie-sv) and bwa (Version: 0.7.17-r1188). The binary file of bwa have been provided in the /ITIPs/scripts/ folder.

```
git clone https://github.com/caixu0518/ITIPs.git
cd ITIPs
chmod u+x *pl 
cd scripts
chmod u+x *
```
# Inputs
Two types of inputs are required for ITIPs
1. Genome fasta. i,e, reference genome, query 1 genome, query 2 genome ......
2. population-scale resequencing reads.  i.e. Sam1_1.fq.gz, Sam1_2.fq.gz ......

# Outputs

Phase 1: the pipeline will generated reference TE insertion and non-reference TE insertion. i.e. reference_TE.insertions.xls, non-reference_TE.insertions.xls
```
Chr     Start   End     Type    SVlen   Upstream        Downstream      Gene    CDS
A10     90093   90797   deletion        704     query1;query2   BraA10g000190.3.1C      BraA10g000200.3.1C      -       -
A10     158873  161879  deletion        3006    query1  BraA10g000350.3.1C      BraA10g000340.3.1C      -       -
A10     161994  162256  deletion        262     query1  BraA10g000350.3.1C      -       -       -
A10     248968  253788  deletion        4820    query2  -       BraA10g000500.3.1C      -       -
A10     252712  253372  deletion        660     query1  -       BraA10g000500.3.1C      -       -
A10     253389  254595  deletion        1206    query1  -       BraA10g000500.3.1C      -       -
A10     253794  254442  deletion        648     query2  -       BraA10g000500.3.1C      -       -
A10     325403  326892  deletion        1489    query1  BraA10g000690.3.1C      -       -       -
A10     325405  325631  deletion        226     query2  BraA10g000680.3.1C      -       -       -
A10     325635  326898  deletion        1263    query2  BraA10g000690.3.1C      -       -       -
A10     329657  332373  deletion        2716    query2  -       BraA10g000690.3.1C      -       -
```

Phase 2: Genotypes of TE insertions in each resequencing genome. i.e. Sam1.refereceTEinsertion and Sam1.Non-refereceTEinsertion.
```
TEindex AB      BC      AC      L       R       Genotype
Dref100 0       0       0       5       0       GG
Dref103 0       0       0       0       0       NA
Dref104 0       0       0       5       0       GG
Dref113 0       0       0       1       0       GG
Dref116 0       0       0       0       0       NA
Dref12  0       0       0       1       0       GG
Dref122 0       0       0       0       0       NA
Dref123 0       0       0       1       0       GG
Dref125 0       0       0       1       0       GG
```
CC indicates that the genotype in the corresponding accession was consistent with the reference genome, and GG indicates that the genotype in the accession was different from the reference genome, NA represents missing loci.


# Usage
There are three main sequential steps to identify and genotype TE insertions, corresponding to 01.Reference_Nonreference_TEinsertion.pl, 02.get_TE_insertions_and_flankingSequences.pl, and 03.TE_insertions_genotype.pl.

Step 1: Identification of reference and non-reference TE insertions between different genomes.

```
perl 01.Reference_Nonreference_TEinsertion.pl  -h

Usage: perl 01.Reference_Nonreference_TEinsertion.pl  -query <query.info.lst>  -ref <reference.info.lst>  -TElib <EDTA.TElib.fa>  -bin <the path to smartie-sv>  -script <the path to scripts>

-query	[required] the query id and query genome files. Two columns (queryName queryGenomeFile).
-ref    [required] the reference information. Three columns (referenceName ReferenceGenomeFile ReferenceGff3)
-TElib  [required] the species TE library
-bin   	[required] the path to smartie-sv. i.e. /10t/caix/src/smartie-sv/bin
-script [required] the path to perl scripts
```

Step 2: extract flanking seuqences of each TE insertion
```
perl 02.get_TE_insertions_and_flankingSequences.pl  -h

Usage: perl 02.get_TE_insertions_and_flankingSequences.pl  -refGenome <ref.fa>  -refName  <ref>   -script <the path to scripts>  reference_TE.insertions.xls  non-reference_TE.insertions.xls

-refGenome	[required] the reference genome in fasta foramt.
-refName	  [required] the reference genome name, same as provided in the first step.
-script   	[required] the path to perl scripts
```

Step 3: Genotype TE insertions using short reads
```
perl 03.TE_insertions_genotype.pl -h

perl 03.TE_insertions_genotype.pl   -Fasta  <ref.referenceTEinsertions_and_flanking1kb.fasta>  -leftRead <Sam1_1.fq.gz>  -rightRead <Sam1_2.fq.gz> -samId <Sam1>  -output <Sam1.refereceTEinsertion>  -script <the path to scripts>  -threads <threads>

-Fasta	[required] the TE insertion and flanking sequences in fasta format
-leftRead	[required] left read
-rightRead	[required] right read
-samId	[required] sample name i.e. Sam1
-output	[required] TE genotype results
-script	[required] the path to scripts
-threads	[optional] threads  default: 6 cores
```

# An example
we recommend that users modified their file format with those we provided in the testData.

```
cd  testData
sh  runMe.sh
cat runMe.sh

#!/bin/bash
path=`pwd`
script="${path}/../scripts/"

perl ../01.Reference_Nonreference_TEinsertion.pl  -query   query.info.lst  \
                                                  -ref     reference.info.lst  \
                                                  -TElib   EDTA.TElib.fa  \
                                                  -bin     /10t/caix/src/smartie-sv/bin  \
                                                  -script  ${script} 

perl ../02.get_TE_insertions_and_flankingSequences.pl  -refGenome ref.fa  \
                                                       -refName   ref  \
                                                       -script    ${script}

perl ../03.TE_insertions_genotype.pl      -Fasta  ref.Non-referenceTEinsertions_and_flanking1kb.fasta  \
                                          -leftRead Sam1_1.fq.gz  \
                                          -rightRead Sam1_2.fq.gz  \
                                          -samId  Sam1 \
                                          -output Sam1.Non-refereceTEinsertion  \
                                          -script  ${script}  \
                                          -threads 10 

perl ../03.TE_insertions_genotype.pl      -Fasta  ref.referenceTEinsertions_and_flanking1kb.fasta  \
                                          -leftRead Sam1_1.fq.gz  \
                                          -rightRead Sam1_2.fq.gz  \
                                          -samId  Sam1 \
                                          -output Sam1.refereceTEinsertion  \
                                          -script  ${script}  \
                                          -threads 10
```
Final results:  
Phase 1: reference_TE.insertions.xls and non-reference_TE.insertions.xls  
Pahse 2: Sam1.refereceTEinsertion and Sam1.Non-refereceTEinsertion

# Citations
