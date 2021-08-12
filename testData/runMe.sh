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
		



